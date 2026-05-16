# PadeTaylor.jl — Public API audit (2026-05-16)

Bead: `padetaylor-qur`. Pre-v1.0 holistic review of the ~30 exported
symbols across 17 modules in `src/`. **Advisory document** — produces
a normalisation recommendation list; actual fixes land in spawned
follow-up beads. No source code was modified during this audit.

---

## 1. Executive summary

The five inconsistencies most worth normalising before v1.0:

| # | Finding | Priority |
|---|---|---|
| 1 | **Step-size kwarg name is `h` in 5 places, `h_max` in 1, `h_path` in 1.** Three names for the same conceptual parameter is the single most user-visible inconsistency in the public surface (`solve_pade(; h_max)` vs `path_network_solve(; h)` vs `lattice_dispatch_solve(; h_path)`). | **v1.0-blocking** |
| 2 | **Newton-iteration cap kwarg name varies: `maxiter` (BVP), `bvp_maxiter` (Dispatcher), `max_iter` (EdgeGatedSolve), `max_steps` (Problems), `max_rescales` (PadeStepper / PathNetwork), `max_steps_per_target` (PathNetwork).** Six different names for "give up after N iterations" across six public functions. | **v1.0-blocking** |
| 3 | **`Diagnostics.quality_diagnose` has no positive-keyword validator** — the docstring promises `sheet`, `tol_well`, `tol_bad`, `n_worst` kwargs, but the generic function declaration `function quality_diagnose end` accepts any kwargs silently. Combined with the extension-based dispatch this means typos go through to a `MethodError` from `DelaunayTriangulation` rather than a clean Rule-1 fail-loud at the package boundary. | cosmetic |
| 4 | **`error("...")` used instead of typed exceptions in 4 sites** (Problems.jl:178, 197; RobustPade.jl:468, 480), one of which is the only fail mode of `solve_pade`'s 1st-order branch (a documented v1 limitation). All other error sites in the project use `throw(ArgumentError|ErrorException|DomainError)`. | cosmetic |
| 5 | **`EdgeReport` / `DiagnosticReport` / `*Solution` naming split** — every other solution-shaped output struct ends in `Solution` (`BVPSolution`, `LatticeSolution`, `DispatcherSolution`, `EdgeGatedSolution`, `PadeTaylorSolution`, `PainleveSolution`, `PathNetworkSolution`, `IVPBVPSolution`). The two diagnostic containers (`DiagnosticReport`, `EdgeReport`) deliberately break the pattern but the naming gap will surprise users reading the export list. | defer |

A non-issue worth flagging in the other direction: **error-message
style discipline is excellent** — 41/56 throws (~73 %) carry an
explicit `Suggestion:` line per Rule 1, and 12 carry an explicit
`detail:` block. The remaining 27 % are uniformly the "validate scalar
range" style (`tol must be positive (got $tol)`) where suggestion is
self-evident.

---

## 2. Audit scope

Files read in full:
- `src/PadeTaylor.jl` (main module, export list — 182 LOC)
- `src/Problems.jl` (261 LOC)
- `src/PathNetwork.jl` (1108 LOC)
- `src/BVP.jl` (592 LOC)
- `src/Dispatcher.jl` (392 LOC)
- `src/EdgeGatedSolve.jl` (353 LOC)
- `src/LatticeDispatcher.jl` (464 LOC)
- `src/RobustPade.jl` (500 LOC)
- `src/Painleve.jl` (375 LOC)
- `src/PainleveSolution.jl` (278 LOC)
- `src/PainleveNamed.jl` (147 LOC)
- `src/PainleveClosedForm.jl` (322 LOC)
- `src/IVPBVPHybrid.jl` (792 LOC)
- `src/Diagnostics.jl` (208 LOC)
- `src/PoleField.jl` (251 LOC)
- `src/EdgeDetector.jl` (262 LOC)
- `src/CoordTransforms.jl` (206 LOC)
- `src/SheetTracker.jl` (322 LOC)
- `src/BranchTracker.jl` (199 LOC)
- `src/PadeStepper.jl` (601 LOC)
- `src/Coefficients.jl` (205 LOC)
- `src/StepControl.jl` (264 LOC)
- `src/LinAlg.jl` (111 LOC)

Skimmed via grep:
- `docs/adr/` (17 ADRs, titles only)

**Exported symbols audited: 39** (re-counted from `src/PadeTaylor.jl:159-180`
including the `painleveplot` extension stub and `PadeTaylorAlg`).

Notable omission: `pii_rational`, `pii_airy`, `piv_entire` are
*defined and exported* from `Painleve` (`src/Painleve.jl:93`) but
**not** re-exported from the main `PadeTaylor` module — the
`PadeTaylor.jl:159-180` export block omits them. See finding (e).3.

---

## 3. (a) Kwarg naming consistency

### Findings

**(a).1 — Step-size kwarg: three names for the same concept.**

The "step length the Padé-Taylor stepper takes" appears under three
distinct kwarg names across public solvers:

| Function | Kwarg name | Default | File:line |
|---|---|---|---|
| `solve_pade` | `h_max` | required | `src/Problems.jl:173` |
| `path_network_solve` | `h` | `0.5` | `src/PathNetwork.jl:377` |
| `dispatch_solve` | `h` | `0.5` | `src/Dispatcher.jl:218` |
| `edge_gated_pole_field_solve` | `h` | `0.5` | `src/EdgeGatedSolve.jl:276` |
| `lattice_dispatch_solve` | `h_path` | `0.5` | `src/LatticeDispatcher.jl:261` |
| `bvp_solve` | (no `h`; uses N collocation nodes) | — | `src/BVP.jl:239` |
| `PadeTaylorAlg` | `h_max` | required | `src/PadeTaylor.jl:140` |

The internal `pade_step!` / `pade_step_with_pade!` / `adaptive_pade_step!`
all use `h` positionally (`src/PadeStepper.jl:267, 300, 553`), so the
public-API split is purely a naming choice, not a semantic difference.

The naming asymmetry is **load-bearing** in `lattice_dispatch_solve`'s
docstring (`src/LatticeDispatcher.jl:240-242`):

> `h_path` — IVP step size for the path-network walker. Forwarded as
> `h_path` to `path_network_solve` (manual-fallback path) **and as
> `h` to `edge_gated_pole_field_solve` (the kwarg name on edge-gated
> is `h`, not `h_path`)**.

That parenthetical is exactly the kind of API smell the audit is
catching. The maintainer's docstring even apologises for it.

**Recommendation:** canonicalise to `h` everywhere. `solve_pade`'s
`h_max` is half-justified (it's a *fixed upper bound* on the
fixed-step driver, not the per-step adaptive `h`), but the name is
still misleading — `solve_pade` does no `min(h_max, h_adapt)`
selection (`src/Problems.jl:202`), it just uses `h_max` as the step
length when `state.z + h_max ≤ z_end`. Rename to `h` and add a
docstring note that `solve_pade` is fixed-step.

**Priority: v1.0-blocking.** This is the single most user-visible
inconsistency. Deprecating with a `Base.@deprecate_binding`-style
shim is the safe path.

---

**(a).2 — Max-iteration kwarg: six different names.**

| Function | Kwarg | Default | File:line |
|---|---|---|---|
| `solve_pade` | `max_steps` | `100_000` | `src/Problems.jl:174` |
| `PadeTaylorAlg` | `max_steps` | `100_000` | `src/PadeTaylor.jl:140` |
| `bvp_solve` (2-arg + 3-arg) | `maxiter` | `10` | `src/BVP.jl:242, 373` |
| `dispatch_solve` | `bvp_maxiter` | `10` | `src/Dispatcher.jl:222` |
| `path_network_solve` | `max_steps_per_target` | `1000` | `src/PathNetwork.jl:394` |
| `path_network_solve` | `max_rescales` | `50` | `src/PathNetwork.jl:384` |
| `adaptive_pade_step!` | `max_rescales` | `50` | `src/PadeStepper.jl:557` |
| `edge_gated_pole_field_solve` | `max_iter` | `500` | `src/EdgeGatedSolve.jl:282` |

These are not all *the same* concept — `max_steps` is "Padé-Taylor
steps until z_end", `max_steps_per_target` is "wedge steps until the
target", `max_rescales` is "FFW adaptive rescales per step",
`maxiter` is "Newton iterations on the BVP slice", `max_iter` is
"region-growth passes" — but the naming **inconsistency is real**:
`maxiter` vs `max_iter` vs `max_steps`-style is purely stylistic.

**Recommendation:** standardise on `max_<noun>` with snake-cased noun
(`max_steps`, `max_iter`, `max_rescales`, `max_steps_per_target`). Rename
`maxiter` → `max_iter` in `bvp_solve` and forward through
`dispatch_solve` (rename `bvp_maxiter` → `bvp_max_iter` for consistency).

**Priority: v1.0-blocking** (BVP has the single un-snake-cased kwarg
in the codebase).

---

**(a).3 — Tolerance kwargs: 5 different names, mostly justified.**

| Function | Kwarg | Default | File:line |
|---|---|---|---|
| `bvp_solve` | `tol` | `eps(T)^(3/4)` | `src/BVP.jl:241` |
| `dispatch_solve` | `bvp_tol` | `nothing` | `src/Dispatcher.jl:221` |
| `dispatch_solve` | `derivative_match_tol` | `1e-7` | `src/Dispatcher.jl:216` |
| `lattice_dispatch_solve` | `bvp_tol` | `nothing` | `src/LatticeDispatcher.jl:265` |
| `path_network_solve` | `adaptive_tol` | `1e-12` | `src/PathNetwork.jl:382` |
| `adaptive_pade_step!` | `adaptive_tol` | `1e-12` | `src/PadeStepper.jl:555` |
| `solve_pole_free_hybrid` | `glue_tol` | `1e-8` | `src/IVPBVPHybrid.jl:442` |
| `robust_pade` | `tol` | `default_tol(T)` | `src/RobustPade.jl:357` |

These names ARE semantically distinct (Newton step-norm, junction
mismatch, adaptive truncation, glue continuity, SV threshold) and the
`<role>_tol` naming pattern is consistent.

**Finding:** `dispatch_solve` and `lattice_dispatch_solve` already use
`bvp_tol` to forward to `bvp_solve`'s `tol`; that pattern works.
However `path_network_solve` exposes `adaptive_tol` directly (not
`step_tol` or `ffw_tol`), which is OK but contributes to the surface
size.

**Recommendation:** leave as-is; the distinct semantics justify
distinct names. Document each in a top-level "kwarg glossary" table
in the package README before v1.0.

**Priority: defer.**

---

**(a).4 — Order kwarg: uniform.**

| Function | Kwarg | Default | File:line |
|---|---|---|---|
| `PadeTaylorProblem` | `order` | `30` | `src/Problems.jl:115` |
| `solve_pade` (via `prob.order`) | — | inherited | `src/Problems.jl` |
| `path_network_solve` | `order` | `prob.order` | `src/PathNetwork.jl:378` |
| `dispatch_solve` | `order` | `prob_ivp.order` | `src/Dispatcher.jl:219` |
| `lattice_dispatch_solve` | `order` | `prob.order` | `src/LatticeDispatcher.jl:262` |
| `edge_gated_pole_field_solve` | `order` | `prob.order` | `src/EdgeGatedSolve.jl:277` |
| `taylor_coefficients_1st/2nd` | `order` (positional) | required | `src/Coefficients.jl:114, 171` |
| `PainleveProblem` | `order` | `30` | per equation, e.g. `src/Painleve.jl:217` |
| `tritronquee`, `hastings_mcleod`, `pii_rational`, `pii_airy`, `piv_entire` | `order` | `30` | `src/PainleveNamed.jl:93, 133`; `src/PainleveClosedForm.jl:206, 267, 302` |
| `adaptive_pade_step!` | `order` (positional) | required | `src/PadeStepper.jl:553` |

This is **well-handled** — every solver inherits `prob.order` by
default, with `order = 30` as the constructor default, and there is
no inconsistency.

**Recommendation:** none.

---

**(a).5 — Verbose / progress kwargs: only two sites, but inconsistent.**

| Function | Kwarg | Default | File:line |
|---|---|---|---|
| `path_network_solve` | `verbose`, `progress_every` | `false`, `500` | `src/PathNetwork.jl:397-398` |
| `edge_gated_pole_field_solve` | `verbose` | `false` | `src/EdgeGatedSolve.jl:283` |

Only two functions take a `verbose` kwarg, but their semantics
differ: `path_network_solve` emits one line every `progress_every`
inner steps via `@printf`; `edge_gated_pole_field_solve` emits a
`@info` per region-growth pass with no granularity knob.

**Recommendation:** if a `verbose` kwarg is exposed, ensure it routes
through `Logging.@info` consistently. The `@printf(stdout, ...)`
pattern in `path_network_solve` bypasses Julia's logging system —
debatable choice that should be either justified in the docstring
("eager-flushed for long progress runs") or normalised.

**Priority: cosmetic.**

---

**(a).6 — Branch / sheet kwargs after the FFW arc.**

| Function | Kwarg(s) | File:line |
|---|---|---|
| `path_network_solve` | `branch_points`, `branch_cut_angles`, `cross_branch`, `initial_sheet`, `grid_sheet` | `src/PathNetwork.jl:386-391` |
| `eval_at_sheet` | `sheet` (positional) | `src/PathNetwork.jl:1028` |
| `eval_at` | (sheet-blind; no kwarg) | `src/PathNetwork.jl:1086` |
| `extract_poles` | (no sheet handling) | `src/PoleField.jl:204` |
| `quality_diagnose` | `sheet=0` | `src/Diagnostics.jl:172` |
| `accumulate_winding`, `sheet_index` | (sheet primitives, positional `branch`) | `src/SheetTracker.jl:304, 320` |

**Inconsistency:** the public surface has *two* sheet conventions in
play:
1. **Multi-branch tuples** (`Vector{Int}` per branch point) used in
   `path_network_solve`'s `grid_sheet` and `eval_at_sheet`'s `sheet`.
2. **Single-branch integer** (`sheet::Int`) used in
   `quality_diagnose` and assumed in `sheet_index`.

These are reconcilable (`quality_diagnose`'s `sheet = 0` is the
principal-branch shorthand) but the type mismatch between
`Vector{Int}` and `Int` for the same conceptual kwarg is the kind of
thing a v1.0 user will trip on.

**Recommendation:** rename `quality_diagnose`'s `sheet::Int = 0` to
`sheet_index::Int = 0` (matches `SheetTracker.sheet_index`'s
semantics — "a single branch's sheet number") OR accept
`AbstractVector{Int}` analogously to `eval_at_sheet`. The latter is
more honest about the multi-branch capability the rest of the API
exposes.

**Priority: cosmetic.**

---

## 4. (b) Return-type consistency

### Findings

**(b).1 — Solution structs: uniform `*Solution` naming, mostly.**

Eight Solution-shaped exports:

| Type | File:line | Defining module |
|---|---|---|
| `PadeTaylorSolution` | `src/Problems.jl:153` | `Problems` |
| `PathNetworkSolution` | `src/PathNetwork.jl:205` | `PathNetwork` |
| `BVPSolution` | `src/BVP.jl:185` | `BVP` |
| `DispatcherSolution` | `src/Dispatcher.jl:168` | `Dispatcher` |
| `EdgeGatedSolution` | `src/EdgeGatedSolve.jl:130` | `EdgeGatedSolve` |
| `LatticeSolution` | `src/LatticeDispatcher.jl:215` | `LatticeDispatcher` |
| `PainleveSolution` | `src/PainleveSolution.jl:80` | `Painleve` |
| `IVPBVPSolution` | `src/IVPBVPHybrid.jl:385` | `IVPBVPHybrid` |

Two diagnostic-shaped exports (NOT `*Solution`):

| Type | File:line | Defining module |
|---|---|---|
| `DiagnosticReport` | `src/Diagnostics.jl:152` | `Diagnostics` |
| `EdgeReport` | `src/Diagnostics.jl:122` | `Diagnostics` |

The `Report` suffix is *deliberate* — these are inspection containers
attached to a `PathNetworkSolution.diagnostics` field, not standalone
solve outputs. The break in pattern is intentional, but the **export
list** (`src/PadeTaylor.jl:163`) places them next to the `*Solution`
exports without comment:

```julia
export path_network_solve, PathNetworkSolution, eval_at, eval_at_sheet
export DiagnosticReport, EdgeReport, quality_diagnose
```

**Recommendation:** keep the `Report` suffix (it correctly signals
"diagnostic, not solve output") but add a short export-list comment
making the distinction explicit. Alternatively rename to
`PathNetworkDiagnostic` and `PathNetworkEdgeDiagnostic` for full
symmetry — but that loses concision.

**Priority: defer** (cosmetic; the distinction is correct).

---

**(b).2 — Solution accessor patterns vary.**

| Solution type | Grid accessor | Callable? | File:line |
|---|---|---|---|
| `PadeTaylorSolution` | `sol.z`, `sol.y` (Vector{Tuple}) | YES, returns `(u, u')` | `src/Problems.jl:217` |
| `PathNetworkSolution` | `sol.grid_z, sol.grid_u, sol.grid_up` | NO (use `eval_at`) | `src/PathNetwork.jl:205` |
| `BVPSolution` | `sol.nodes_z, sol.u_nodes` | YES, returns `(u, u')` | `src/BVP.jl:486` |
| `DispatcherSolution` | `sol.grid_z, sol.grid_u, sol.grid_up` | NO | `src/Dispatcher.jl:168` |
| `EdgeGatedSolution` | `sol.u_grid` (Matrix), `sol.xs, sol.ys` | NO | `src/EdgeGatedSolve.jl:130` |
| `LatticeSolution` | `sol.grid_z, sol.grid_u, sol.grid_up` (Matrix) | NO | `src/LatticeDispatcher.jl:215` |
| `PainleveSolution` | `grid_values(sol)` accessor | YES *if* raw is callable | `src/PainleveSolution.jl:132` |
| `IVPBVPSolution` | (no grid; uses underlying PFS grids) | YES, in BVP sector | `src/IVPBVPHybrid.jl:732` |

**Findings:**
- The "grid_z / grid_u / grid_up" idiom is consistent across 4 types
  (PathNetwork, Dispatcher, Lattice, EdgeGated *partially*).
- `EdgeGatedSolution` uses `u_grid` (Matrix) instead of `grid_u`
  (Vector), reflecting the 2D-lattice nature but breaking the naming
  pattern.
- `LatticeSolution` uses `grid_u` as `Matrix{Complex{T}}` (2D),
  whereas other `grid_u`s are `Vector{Complex{T}}` (1D scatter).
  Same name, different shape.
- Direct field access vs functional accessor: `PainleveSolution`
  consciously uses `grid_values(sol)` to abstract the underlying raw
  type. Other Solution structs expose fields directly.

**Recommendation:** before v1.0, document the access pattern in each
Solution struct's docstring with a `## Accessing the output` section.
For `LatticeSolution.grid_u` (Matrix) vs `PathNetworkSolution.grid_u`
(Vector), keep the name but **highlight** the shape difference in
both docstrings.

**Priority: cosmetic.**

---

**(b).3 — Top-level solvers all return Solution structs.**

Verified: every top-level public solver returns a `*Solution` struct,
never a raw tuple / Vector / Dict. The only "raw" return is
`extract_poles` returning `Vector{<:Complex}` (correct — it is a
single-purpose post-hoc query, not a solver).

`eval_at` and `eval_at_sheet` return `Tuple{Complex, Complex}` —
correct for per-point query.

**Recommendation:** none.

---

**(b).4 — Predicate / mutator naming.**

The exported surface has:

- **One Bool-returning predicate**: none exported, but several
  internal: `any_cut_crossed` (`src/BranchTracker.jl:130`),
  `segment_crosses_cut` (`src/BranchTracker.jl:112`) — both internal
  to `BranchTracker`, NOT exported. Good — Bool predicates as public
  API would warrant an `is_*` prefix discussion, but none exist.

- **Symbol-returning classifiers**: none exported. The `region_tag`
  values inside `LatticeSolution` (`:ivp`, `:bvp`, `:bvp_fail`,
  `:ivp_only`) are field values, not function returns.

- **Mutator (`!`-suffix) functions**: only three, all in
  `PadeStepper`: `pade_step!`, `pade_step_with_pade!`,
  `adaptive_pade_step!`. None are exported from the main module.

**Recommendation:** the `!`-suffix internal-only convention is
correct (mutators are stepping primitives, not user-facing). No
predicate-naming issue exists because no Bool predicates are
exported. Leave as-is.

**Priority: defer.**

---

## 5. (c) Error message style

### Findings

**(c).1 — Sample sites and pattern conformance.**

Sampling 10 throw sites across modules:

| File:line | Function | Has `Suggestion:` | Has `detail:` | Has `(got …)` |
|---|---|---|---|---|
| `src/Problems.jl:116-118` | `PadeTaylorProblem` | YES | NO | YES |
| `src/Problems.jl:175-177` | `solve_pade` | YES | NO | YES |
| `src/BVP.jl:246-249` | `bvp_solve` | YES | NO | YES |
| `src/BVP.jl:250-251` | `bvp_solve` | NO | NO | YES |
| `src/BVP.jl:262-263` | `bvp_solve` | NO | NO | YES |
| `src/BVP.jl:332-338` | `bvp_solve` | YES (multi) | NO | NO |
| `src/PathNetwork.jl:417-420` | `path_network_solve` | NO | NO | YES |
| `src/PathNetwork.jl:526-529` | `path_network_solve` | NO | NO | YES |
| `src/PathNetwork.jl:586-590` | `path_network_solve` | YES | NO | NO |
| `src/PathNetwork.jl:870-873` | `path_network_solve` (R helper) | NO | NO | YES |
| `src/EdgeGatedSolve.jl:285-287` | `edge_gated_pole_field_solve` | NO | YES | YES |
| `src/EdgeGatedSolve.jl:288-289` | `edge_gated_pole_field_solve` | NO | NO | NO |
| `src/LatticeDispatcher.jl:270-272` | `lattice_dispatch_solve` | NO | YES | YES |
| `src/LatticeDispatcher.jl:280-282` | `lattice_dispatch_solve` | NO | NO | YES |
| `src/Painleve.jl:164-168` | `PainleveProblem` | YES | NO | NO |
| `src/Painleve.jl:248-251` | `_build_transformed` | YES | NO | NO |

**Pattern observation:** `Suggestion:` is used 41 times (74% of the
123 throw-with-message sites), `detail:` is used 12 times (10%).
There are at least three "styles":
- The full **prose** style with explicit `Suggestion:` (most common).
- The terse **range-validate** style (`tol must be positive (got $tol)`)
  with no suggestion (~20 sites; e.g. `src/BVP.jl:251`,
  `src/PathNetwork.jl:430`).
- The **mixed-prose-detail** style using `detail:` instead of
  `Suggestion:` (10 sites, mostly in EdgeDetector / EdgeGatedSolve /
  LatticeDispatcher — these were authored later in the project).

**(c).2 — `error(...)` calls without an exception type.**

Found 4 sites:

| File:line | Function | Context |
|---|---|---|
| `src/Problems.jl:178-182` | `solve_pade` | "1st-order branch not implemented in v1" |
| `src/Problems.jl:197-201` | `solve_pade` | "did not reach z_end after max_steps" |
| `src/RobustPade.jl:468-470` | `_trim_and_normalise` | "all denominator coefficients below tol" |
| `src/RobustPade.jl:480` | `_trim_and_normalise` | "post-trim b is empty" |

All others use `throw(ArgumentError|ErrorException|DomainError)`. The
4 bare `error()` calls produce `ErrorException` at runtime so the
type signature is fine, but stylistically inconsistent.

**Recommendation:** replace bare `error(...)` with explicit
`throw(ErrorException(...))` for all four sites. Mechanically trivial;
the only judgement call is `Problems.jl:178` ("1st-order branch not
implemented in v1") — this is arguably a `MethodError` (the method
isn't defined) but `ErrorException` with a `Suggestion:` line is the
honest fail-loud per Rule 1.

**(c).3 — `(got X)` formatting consistency.**

Most sites use `(got $value)` placed immediately after the validation
clause: `"order must be ≥ 2 (got $order)"`. A handful use embedded
formatting:
- `src/EdgeGatedSolve.jl:286-287`: `"need length(xs), length(ys) ≥ 3 (got $nx, $ny)"`.
- `src/PainleveClosedForm.jl:99`: `"n must be 1, 2, or 3 (got $n)."`.

These are consistent within local style; no normalisation needed.

**Recommendation overall for (c):** convert the 4 bare `error()`
calls to typed `throw(ErrorException(...))`. Optionally normalise the
terse range-validate sites (~20) to include a `Suggestion:` line on
the next pass — but this is cosmetic.

**Priority: cosmetic, but the bare `error()`-replacement is a
mechanical pre-v1.0 cleanup worth doing.**

---

## 6. (d) Bang-suffix usage

### Findings

**(d).1 — `!`-suffix functions in `src/`.**

Only three:

| Function | File:line | Mutates first arg? |
|---|---|---|
| `pade_step!` | `src/PadeStepper.jl:267` | YES — mutates `state::PadeStepperState{T}` |
| `pade_step_with_pade!` | `src/PadeStepper.jl:300` | YES — mutates `state::PadeStepperState{T}` |
| `adaptive_pade_step!` | `src/PadeStepper.jl:553` | YES — mutates `state::PadeStepperState{T}` |

All three correctly mutate their first argument and return it (plus
the local Padé / metadata). **None are exported** from the main
`PadeTaylor` module (only `PadeStepperState` would need to be too —
neither is in `src/PadeTaylor.jl:159-180`).

**(d).2 — Non-`!` functions that mutate.**

Spot-checked driver functions for mutation:
- `solve_pade` mutates a private `state` local, returns a fresh
  `PadeTaylorSolution` — correct.
- `path_network_solve` allocates fresh arrays — correct.
- `bvp_solve` mutates a private `u_int` local — correct.

No exported non-`!` function silently mutates a caller's argument.

**Recommendation:** none — the convention is honoured cleanly.

**Priority: defer** (already correct).

---

## 7. (e) Capitalization / case style

### Findings

**(e).1 — Function names: snake_case throughout.**

All exported function names conform: `path_network_solve`,
`extract_poles`, `eval_at`, `eval_at_sheet`, `lattice_dispatch_solve`,
`taylor_coefficients_1st`, etc. No camelCase or PascalCase outliers.

**(e).2 — Type names: PascalCase throughout.**

All exported type names conform: `PadeTaylorProblem`,
`PadeTaylorSolution`, `PathNetworkSolution`, `BVPSolution`,
`DiagnosticReport`, etc.

One curiosity: `bvp_solve` is snake_case but `BVPSolution` is
all-caps acronym + PascalCase. This is the Julia-idiomatic choice for
acronyms (matches `LAPACK`, `FFT`, `SVD` precedents). Fine.

**(e).3 — Export-list omission.**

The `Painleve` module declares `pii_rational`, `pii_airy`,
`piv_entire` in its export list (`src/Painleve.jl:93`) but they are
NOT re-exported from the top-level `PadeTaylor` module
(`src/PadeTaylor.jl:175-178`). Per CLAUDE.md Rule 1, this is silent
behaviour — calling `using PadeTaylor; pii_rational(1)` would fail
with `UndefVarError` despite the `Painleve` constructor existing in
the package.

The closest analog (`tritronquee`, `hastings_mcleod`) ARE re-exported
on `src/PadeTaylor.jl:177`.

**Recommendation:** add `pii_rational, pii_airy, piv_entire` to the
top-level export list. Single-line fix.

**Priority: v1.0-blocking** (otherwise these constructors are
discoverable only via `PadeTaylor.Painleve.pii_rational`, which
violates the "named-transcendent constructors are the discoverable
home" promise made in `PainleveClosedForm.jl:1-13`).

---

## 8. Additional checks

### 8.1 — Symbol-valued kwargs

**Findings.** Symbol kwargs are used in:

| Function | Kwarg | Allowed values | File:line |
|---|---|---|---|
| `path_network_solve` | `step_selection` | `:min_u`, `:steepest_descent` | `src/PathNetwork.jl:380` |
| `path_network_solve` | `step_size_policy` | `:fixed`, `:adaptive_ffw` | `src/PathNetwork.jl:381` |
| `pole_field_mask` | `level` | `:auto` OR a numeric | `src/EdgeDetector.jl:216` |
| `edge_gated_pole_field_solve` | `edge_level` | `:auto` OR a numeric | `src/EdgeGatedSolve.jl:278` |
| `lattice_dispatch_solve` | `edge_level` | `:auto` OR a numeric | `src/LatticeDispatcher.jl:263` |
| `PainleveProblem` | (positional `equation`) | `:I`, `:II`, `:III`, `:IV`, `:V`, `:VI` | `src/Painleve.jl:202` |
| `PainleveProblem(:VI; frame)` | `frame` | `:transformed`, `:transformed_eta` | `src/Painleve.jl:277` |
| `tritronquee` | (positional `equation`) | `:I` only | `src/PainleveNamed.jl:92` |
| `hastings_mcleod` | `branch` | `:positive`, `:negative` | `src/PainleveNamed.jl:132` |
| `piv_entire` | (positional `kind`) | `:minus_2z`, `:minus_two_thirds_z` | `src/PainleveClosedForm.jl:302` |
| `LatticeSolution.region_tag` | (field values) | `:ivp`, `:bvp`, `:bvp_fail`, `:ivp_only` | `src/LatticeDispatcher.jl:370` |
| `PainleveProblem.frame` | (struct field) | `:direct`, `:transformed`, `:transformed_eta` | `src/Painleve.jl:147` |

**Findings:**
1. The `edge_level` vs `level` split (same concept, two names — `level`
   in `pole_field_mask`, `edge_level` in the higher-level drivers) is
   a minor wart. (`src/EdgeDetector.jl:216` vs
   `src/EdgeGatedSolve.jl:278`.)
2. Every Symbol-accepting function has an explicit
   `<value> ∈ (...) || throw(ArgumentError(...))` validator with
   `(got :$sym)` formatting. Excellent.
3. Symbol values are **not namespaced** (no `:PadeTaylor.fixed` etc).
   That is the Julia convention — accept the bare keyword.
4. **All documented.** Every Symbol value appears in the calling
   function's docstring with semantics, except `:transformed_eta`
   (`src/Painleve.jl:277-294` — documented in line, missing from the
   PainleveProblem top-level signature docstring at line 188-201).

**Recommendation:**
- Rename `EdgeDetector.pole_field_mask`'s `level` → `edge_level` for
  consistency with downstream consumers. (Or rename `edge_level` →
  `level` everywhere; the former is closer to FW2011 paper terminology
  per `EdgeDetector` module docstring.)
- Add `:transformed_eta` to the `PainleveProblem(:VI; ...)` keyword
  signature in the top-level docstring.

**Priority: cosmetic.**

---

### 8.2 — Required vs optional positional args (4+ smell)

**Findings.** Functions with 4+ positional args:

- `bvp_solve(f, ∂f_∂u, z_a, z_b, u_a, u_b; ...)` — **6 positional** (`src/BVP.jl:239`)
- `bvp_solve(f, ∂f_∂u, ∂f_∂up, z_a, z_b, u_a, u_b; ...)` — **7 positional** (`src/BVP.jl:370`)
- `dispatch_solve(prob_ivp, f_bvp, ∂f_∂u_bvp, segments; ...)` — **4 positional** (`src/Dispatcher.jl:212`)
- `lattice_dispatch_solve(prob, bvp_f, bvp_∂f_∂u, xs, ys; ...)` — **5 positional** (`src/LatticeDispatcher.jl:257`)

The BVP signature is the worst: 7 positional args (`f, ∂f_∂u, ∂f_∂up,
z_a, z_b, u_a, u_b`). The function-callable trio + the four BVP
endpoints make positional sense (they are all required, and the
function-callables go first by convention), but ordering is
memory-load. The maintainer's docstring acknowledges this by spelling
out each arg's role inline.

**Recommendation:** for v1.0, consider a `BVPProblem` constructor
mirroring `PadeTaylorProblem` to bundle the 4 BVP endpoints into a
single struct, leaving `bvp_solve(prob_bvp, f, ∂f_∂u[, ∂f_∂up]; ...)`
as the public driver. This matches the `solve_pade(prob)` pattern
already shipped.

**Priority: cosmetic, but a clean pre-v1.0 refactor opportunity.**

---

### 8.3 — Constructor-vs-solver patterns across Painlevé family

**Findings.** Three patterns coexist:

1. **Equation-symbol keyword constructor**: `PainleveProblem(:I; u0,
   up0, zspan, order)` (`src/Painleve.jl:202`). Required-kwarg
   validator throws on missing/extra (`src/Painleve.jl:160-176`).
2. **Named-IC constructor**: `tritronquee(:I; zspan, order)`
   (`src/PainleveNamed.jl:92`); `hastings_mcleod(; branch, zspan,
   order)` (`src/PainleveNamed.jl:132`). Encapsulate a literature-
   pinned IC.
3. **Closed-form family constructor**: `pii_rational(n; zspan, order)`
   (`src/PainleveClosedForm.jl:206`); `pii_airy(n; θ, zspan, order)`
   (`src/PainleveClosedForm.jl:266`); `piv_entire(kind; zspan, order)`
   (`src/PainleveClosedForm.jl:302`). Derive IC from closed form.

All three return `PainleveProblem` (the `name::Union{Symbol,Nothing}`
field tags which family) — the unifying type design is sound. The
**naming is inconsistent**:
- `tritronquee` (lowercase noun)
- `hastings_mcleod` (lowercase, underscored)
- `pii_rational` (lowercase, prefixed by equation tag)
- `pii_airy` (same)
- `piv_entire` (same)

The first two are *named transcendents* (one specific solution,
named after authors); the last three are *families* (parametrised
classes). The prefix on the family constructors (`pii_`, `piv_`) is a
disambiguator that the named-transcendent constructors lack —
arguably `pii_rational` should be `painleveII_rational` for full
consistency, but the `pii_` shorthand is a reasonable compromise.

**Recommendation:**
1. **Fix the export omission** (see 7.e.3) so `pii_rational` etc are
   discoverable.
2. Document the constructor taxonomy (`PainleveProblem` general; named
   transcendents; closed-form families) in the package README before
   v1.0.

**Priority: v1.0-blocking** (export omission); **defer** (taxonomy
documentation).

---

## 9. Proposed follow-up beads

Orchestrator will file these as separate beads (audit doc is
advisory; this section is a recommendation table).

| # | Title | Scope | Priority |
|---|---|---|---|
| 1 | Canonicalise step-size kwarg name to `h` across all public solvers | Rename `h_max` (`solve_pade`, `PadeTaylorAlg`) and `h_path` (`lattice_dispatch_solve`) to `h`; provide `Base.@deprecate_binding` shims for one release; update docs. | v1.0-blocking |
| 2 | Normalise max-iteration kwarg names (`maxiter` → `max_iter`) | Rename `bvp_solve`'s `maxiter` to `max_iter`; cascade through `dispatch_solve`'s `bvp_maxiter` → `bvp_max_iter`. Deprecate-shim the old names. | v1.0-blocking |
| 3 | Export `pii_rational`, `pii_airy`, `piv_entire` from `PadeTaylor` | One-line addition to `src/PadeTaylor.jl` export block. | v1.0-blocking |
| 4 | Replace bare `error(...)` with typed `throw(ErrorException(...))` | Four sites: `src/Problems.jl:178, 197`; `src/RobustPade.jl:468, 480`. Mechanical. | cosmetic |
| 5 | Add `Suggestion:` line to the ~20 terse range-validate throws | Touch every `throw(ArgumentError("X must be Y (got $Z)"))` that lacks a `Suggestion:` line. Use the BVP-style `Suggestion: (a) … (b) …` template where multi-cause. | cosmetic |
| 6 | Introduce `BVPProblem` constructor + thin-arg `bvp_solve` | Mirror the `PadeTaylorProblem` / `solve_pade` shape: `BVPProblem(z_a, z_b, u_a, u_b; ...)` bundles endpoints; `bvp_solve(bvp_prob, f, ∂f_∂u[, ∂f_∂up]; ...)` takes the bundle. Reduces 7-positional to 3-positional. | cosmetic |
| 7 | Rename `EdgeDetector.pole_field_mask`'s `level` → `edge_level` | One-arg-name change; deprecate `level`. Consistent with downstream consumers. | cosmetic |
| 8 | Unify `quality_diagnose(; sheet)` to accept multi-branch sheet tuples | Currently `Int = 0`; align with `eval_at_sheet`'s `AbstractVector{<:Integer}`. Or rename to `sheet_index::Int` for "I am a single-branch wrapper". | cosmetic |
| 9 | Document Symbol-valued kwargs and accepted values in a glossary | New `docs/src/api_glossary.md` section listing every `:value`-accepting kwarg, its accepted set, and its semantics. README link. | defer |
| 10 | Document constructor taxonomy for Painlevé family | README section explaining `PainleveProblem` (general), named transcendents, closed-form families; how `name::Symbol` field distinguishes them. | defer |
| 11 | Audit `LatticeSolution.grid_u` (Matrix) vs `PathNetworkSolution.grid_u` (Vector) shape collision | Document in both docstrings; consider rename `grid_u_matrix` for the 2D case. | defer |
| 12 | Verbose-mode routing through `Logging.@info` consistently | `path_network_solve` uses `@printf(stdout, ...)`; `edge_gated_pole_field_solve` uses `@info`. Pick one. | defer |

**Bead-count summary:** 12 proposed; 3 v1.0-blocking, 5 cosmetic, 4
defer.

---

## 10. Things that are already consistent

(Audit validation — these are positive findings, NOT complaints.)

1. **All solvers return `*Solution` (or named diagnostic-report)
   structs** — no raw tuple/array returns from public drivers.
2. **`!`-suffix discipline is clean** — the three `!` functions all
   mutate their first argument (`PadeStepperState`); no non-`!`
   functions silently mutate exported state.
3. **PascalCase / snake_case discipline is perfect** for exported
   names.
4. **`order` kwarg threading from `PadeTaylorProblem` → downstream
   solvers is uniform** — every consumer defaults to `prob.order`,
   and `PadeTaylorProblem`'s default is `30` (FW 2011 canonical).
5. **Symbol-valued kwargs are uniformly validated** with explicit
   `∈ (allowed) || throw(...)` patterns and `(got :$sym)` formatting.
6. **74 % of throws include explicit `Suggestion:` text** per CLAUDE.md
   Rule 1 — well above any reasonable threshold for "the rule is
   honoured in practice".
7. **`PainleveProblem` keyword-validator** (`src/Painleve.jl:160-176`)
   is the gold standard: throws on extras AND missing required, with
   a `Suggestion:` line pointing to the canonical form. Every
   constructor in the codebase should be checked against this pattern.
8. **Per-segment Padé store discipline** — `PathNetworkSolution`,
   `PadeTaylorSolution`, and the underlying `BVPSolution` all expose
   the per-node Padé approximants, which `PoleField.extract_poles`
   reads back via a single shared `_extract_poles_core` (`src/PoleField.jl:119`).
   Excellent code-reuse design.
9. **Extension-pattern adherence** (ADR-0003) — `painleveplot`,
   `quality_diagnose`, `PadeTaylorAlg` are all empty-generic /
   data-struct declarations in the main module with implementation
   in `ext/`; no algorithmic logic leaks into extensions.
10. **Backward-compat constructors on `PathNetworkSolution`**
    (`src/PathNetwork.jl:221-230`) maintain pre-ADR-0013 / pre-ADR-0016
    test compatibility via the 9-arg and 10-arg variants. Disciplined.

---

## 11. Audit limitations

- The audit did NOT exercise the test suite (`test/runtests.jl`) —
  per the bead instruction "DO NOT run julia".
- Internal modules (`StepControl`, `Coefficients`, `LinAlg`,
  `PadeStepper`, `BranchTracker`) are *not* exported and were
  audited only for kwarg/error patterns affecting public surface
  behaviour, not for completeness.
- The `ext/PadeTaylorDiagnosticsExt.jl` and
  `ext/PadeTaylorCommonSolveExt.jl` were NOT read in this audit
  (they appear in `git status` as untracked / WIP). Conclusions
  about `quality_diagnose` are based on the in-module generic
  declaration only.
- ADR-0017's `strict` kwarg semantics on `lattice_dispatch_solve`
  (fail-soft mode for non-convergent BVP rows) were audited for
  signature consistency only; the implementation logic was not
  cross-checked against the ADR.

---

*End of audit document. Output is advisory; orchestrator should file
the 12 proposed follow-up beads with the priority tags above.*
