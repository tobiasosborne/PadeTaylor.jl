# ADR-0006 — `PainleveProblem` convenience layer

**Status**: Proposed (2026-05-14) | **Bead**: `padetaylor-avl`
**Context**: Design discussion 2026-05-14 — "is PadeTaylor.jl complete
enough to provide a general-purpose Painlevé function?"  The honest
answer (recorded below) is that the *solver* substrate is general, but
per-equation *problem construction* is missing and lopsided, and a
uniform `painleve()` value-function is the wrong abstraction.

## Context

PadeTaylor.jl's solving layer is equation-agnostic: `solve_pade`,
`path_network_solve`, `bvp_solve`, `dispatch_solve`, and
`lattice_dispatch_solve` all take a `PadeTaylorProblem` wrapping an
arbitrary 2nd-order RHS `f(z, u, u')`.  The FW-figure reproduction work
(worklogs 022–025, 8 figures) drove PI through this substrate without
touching `src/`.

But constructing the *problem* for a given Painlevé equation is today
ad-hoc and unevenly covered:

| equation | RHS factory in `src/` | coordinate transform | movable singularities |
|----------|------------------------|----------------------|-----------------------|
| PI   | none — `6u²+z` written inline in tests/figures | not needed | poles |
| PII  | **none** | not needed | poles |
| PIII | `CoordTransforms.pIII_transformed_rhs(α,β,γ,δ)` | exp, shipped | poles + fixed branch point at `z=0` |
| PIV  | **none** | not needed | poles |
| PV   | `CoordTransforms.pV_transformed_rhs(α,β,γ,δ)` | exp, shipped | poles + fixed branch point at `z=0` |
| PVI  | `SheetTracker.pVI_transformed_rhs(α,β,γ,δ)` | exp (reuses PV's), shipped | poles + fixed branch points at `z=0,1,∞` |

So PIII/PV/PVI have factories; PI/PII/PIV do not.  Worse, the three
existing factories are *not* discoverable as "the Painlevé equations" —
they live in transform-specific Tier-4/5 modules with transformed-frame
signatures.  A user who wants "PII with parameter α and these ICs" has
no obvious entry point.

### Why not a uniform `painleve(eq, z, …)` value-function

Rejected.  The Painlevé transcendents are **not** functions one looks
up like `sin`: every equation has *infinitely many* solutions, and the
interesting ones (tronquée, tritronquée, rational, Hastings–McLeod, …)
are singled out by initial/boundary conditions, not by a canonical
default.  A `painleve(:I, z)` returning "the" value would have to
silently pick a solution — a direct violation of CLAUDE.md Rule 1
(fail loud; never paper over a required decision).  The caller *must*
state which solution they want.

## Decision

Add a per-equation **problem builder**, not a value-function.

**File**: `src/Painleve.jl` (≤200 LOC; split to `src/PainleveRHS.jl`
for the three new RHS factories if it exceeds the budget).

**The wrapper struct**:

```julia
struct PainleveProblem{P}
    equation  :: Symbol                 # :I … :VI
    params    :: NamedTuple             # (), (α=…,), (α=…,β=…,γ=…,δ=…)
    problem   :: P                      # the underlying PadeTaylorProblem
    frame     :: Symbol                 # :direct (I,II,IV) | :transformed (III,V,VI)
    to_frame  :: Function               # z-frame → solve-frame  (identity if :direct)
    from_frame:: Function               # solve-frame → z-frame  (identity if :direct)
end
```

**The constructor** (one method per equation, dispatched on a `Val`
of the equation symbol so unknown symbols fail at construction):

```julia
PainleveProblem(:I;   u0, up0, zspan, order = 30)
PainleveProblem(:II;  α, u0, up0, zspan, order = 30)
PainleveProblem(:IV;  α, β, u0, up0, zspan, order = 30)
PainleveProblem(:III; α, β, γ, δ, z0, u0, up0, zspan, order = 30)
PainleveProblem(:V;   α, β, γ, δ, z0, u0, up0, zspan, order = 30)
PainleveProblem(:VI;  α, β, γ, δ, z0, u0, up0, zspan, order = 30)
```

What the constructor does, and only this:

  1. **Selects the RHS.**  For PI/PII/PIV, builds the direct RHS
     closure from the canonical form + parameters.  For PIII/PV/PVI,
     delegates to the *existing* `CoordTransforms` / `SheetTracker`
     factories — `Painleve.jl` wraps them, it does not duplicate them.
  2. **Handles the coordinate frame.**  For PIII/PV/PVI the caller
     supplies ICs and domain in the natural `z`-frame; the constructor
     applies `pIII_z_to_ζ` / `pV_z_to_ζ` (PVI reuses PV's, FFW 2017
     §2.2) to land the IC in the solve-frame and stores both the
     forward and inverse maps so results can be mapped back.  For
     PI/PII/PIV the frame is `:direct` and both maps are `identity`.
  3. **Bookkeeps parameters.**  Stores `equation` + `params` so the
     problem is self-describing (printing, dispatch, provenance).
  4. **Builds the `PadeTaylorProblem`.**  In the solve-frame, with the
     given `order`.

What the constructor does **not** do:

  - It does **not** solve anything — solving stays the job of
    `solve_pade` / `path_network_solve` / `bvp_solve` / the
    dispatchers.
  - It does **not** pick a solution — ICs are required, never
    defaulted (see Rejected Alternatives).
  - It does **not** add automatic Riemann-sheet or branch-cut routing
    for PIII/PV/PVI walks — that remains the deliberately-deferred
    Tier-5 work (`docs/figure_catalogue.md` T5 PARTIAL; `SheetTracker`
    ships sheet *bookkeeping*, not sheet-constrained *routing*).

**Forwarding methods** (the ergonomic payoff).  v1 ships thin methods

```julia
solve_pade(pp::PainleveProblem; kwargs...)
path_network_solve(pp::PainleveProblem, grid; kwargs...)
bvp_solve(pp::PainleveProblem, …; kwargs...)
```

that (a) unwrap to `pp.problem`, (b) for `:transformed` problems map
the input grid/domain `z → solve-frame` on the way in and the solution
`solve-frame → z` on the way out via `pp.to_frame` / `pp.from_frame`,
and (c) **thread `pp.problem.order` explicitly** into
`path_network_solve` — see Constraints.

## The six canonical forms — where ground truth lives

Per CLAUDE.md Law 1, the implementation opens these at coding time; the
ADR only records *which* reference is authoritative for each:

  - **PI** — `references/markdown/FW2011_painleve_methodology_JCP230/…`
    : `u'' = 6u² + z`.  (Already exercised throughout the figure work.)
  - **PII** — `references/markdown/FW2014_second_PII_exploration_FoCM14/…`
    and `…/FW2015_imaginary_PII_PhysicaD309/…` : `u'' = 2u³ + zu + α`.
  - **PIV** — `references/markdown/ReegerFornberg2014_PIV_fundamental_domain_PhysicaD280/…`
    : the form with the `(u')²/(2u)` term — note the `u = 0`
    singularity (see Fail-loud).
  - **PIII, PV** — `CoordTransforms.jl`'s top docstring already cites
    FFW 2017 §2.1 (`FFW2017_…md:39-49`) verbatim; the transformed RHS
    factories are shipped and tested.
  - **PVI** — `SheetTracker.jl`'s docstring cites FFW 2017 §2.2
    (`FFW2017_…md:144`); `pVI_transformed_rhs` is shipped and tested.

The three *new* factories needed are PI, PII, PIV — all direct-frame,
all with only movable poles, so they drop straight into the existing
path-network machinery once written.

## Constraints carried from existing code

  - **`path_network_solve` ignores `prob.order`** (bead
    `padetaylor-9xf`).  It has its own `order::Integer = 30` kwarg and
    never reads the `PadeTaylorProblem`'s `order` field.  The
    `PainleveProblem` forwarding method for `path_network_solve` MUST
    pass `order = pp.problem.order` explicitly, or the carefully
    constructed `order` is silently dropped.  If `padetaylor-9xf` is
    fixed first (kwarg defaults to `prob.order`), this constraint
    relaxes — but the ADR records it so the implementation cannot
    re-step on the rake worklog 025 already stepped on.

  - **`PadeTaylorProblem` rejects `zspan` with equal endpoints.**  The
    `PainleveProblem` constructor must pick a non-degenerate `zspan`
    for the underlying problem even when the caller only cares about
    the IC point (the path-network use-case) — same workaround the
    figure scripts use (`(0.0, 10.0)` placeholder).

## Fail-loud behaviour (CLAUDE.md Rule 1)

  - Unknown equation symbol → `ArgumentError` listing `:I … :VI`.
  - Missing a required parameter for the equation (e.g. `:II` without
    `α`, `:V` without `δ`) → `ArgumentError` naming the equation, the
    missing parameter, and the canonical form's reference.
  - Supplying a parameter the equation does not take (e.g. `α` to
    `:I`) → `ArgumentError`, not silent ignore.
  - **PIV `u = 0`**: the canonical PIV RHS has a `β/u` term and a
    `(u')²/(2u)` term — singular at `u = 0`.  The constructor cannot
    prevent a walk from approaching `u = 0`; it documents the hazard
    in the docstring and the RHS closure throws a `DomainError` with a
    `suggestion` if evaluated at `u = 0`, rather than returning `Inf`.
  - PIII/PV/PVI with a `z`-frame IC at `z = 0` (the branch point) →
    `ArgumentError`: the transform is singular there; suggest an IC
    off the branch point.

## Scope

**v1 (this bead):**
  - `src/Painleve.jl` — the wrapper struct + six constructors.
  - Three new direct-frame RHS factories (PI, PII, PIV).
  - Forwarding methods for `solve_pade`, `path_network_solve`,
    `bvp_solve`.
  - Re-export `PainleveProblem` from the umbrella module.
  - Tests: one constructor-correctness test per equation (the built
    RHS matches a hand-derived value at a sample point), the
    fail-loud guards, and an end-to-end "PII through the path network
    gives a finite pole field" smoke test.  Mutation-prove the
    PI/PII/PIV factory tests (perturb a coefficient → RED).

**Deferred (separate beads, not v1):**
  - Automatic sheet-aware / branch-cut-avoiding routing for PIII/PV/PVI
    walks — the standing Tier-5 gap.
  - Named special-solution constructors (`tritronquee(:I)`,
    `hastings_mcleod(:II)`, …) — a convenience on top of this layer
    once it exists; still requires the caller to *opt in* to a named
    solution, so it does not reintroduce the rejected value-function.
  - Parameter-space surveys (the FW 2014/2015 IC-plane figures).

## Rejected alternatives

  1. **A `painleve(eq, z; …)` value-function** — rejected above:
     no canonical solution, would force a silent choice.
  2. **Make `PainleveProblem <: PadeTaylorProblem` (subtype) so it
     drops into solvers directly** — rejected: `PadeTaylorProblem` is
     a concrete struct, and the wrapper genuinely needs to carry the
     transform pair for III/V/VI.  Forwarding methods are cleaner than
     a type-hierarchy contortion.
  3. **Default ICs per equation** (e.g. default `:I` to the
     tritronquée) — rejected: it buries a mathematical choice in an
     API default.  Named special-solution constructors (deferred) are
     the right place for "I want *that* solution" — explicit, opt-in.
  4. **Put the PI/PII/PIV factories in `CoordTransforms.jl`** —
     rejected: `CoordTransforms` is specifically the *branch-point
     transform* module; PI/PII/PIV need no transform.  A dedicated
     `Painleve.jl` is the honest home and the discoverable entry point.

## Consequences

  - **Positive**: one discoverable entry point for all six equations;
    PI/PII/PIV become first-class instead of inline-RHS folklore;
    III/V/VI's transform bookkeeping stops leaking into user code;
    results auto-map back to the `z`-frame.
  - **Positive**: the layer is thin and additive — no change to the
    solver substrate, so the 1432-test suite is unaffected except for
    the new `Painleve.jl` tests.
  - **Negative / honest**: it does *not* make PadeTaylor.jl a
    "Painlevé library" in the full sense — the hard part (sheet-aware
    routing for III/V/VI) stays deferred.  The ADR is explicit that
    v1 is a *problem-construction* convenience, not a closing of the
    Tier-5 gap.

## References

  - **CLAUDE.md Rule 1** — fail loud; never default a required
    decision.  The core reason this is a builder, not a function.
  - **ADR-0001** — four-layer architecture; `PainleveProblem` sits
    just above the `Problems` driver layer as a problem-construction
    convenience.
  - **`src/CoordTransforms.jl`**, **`src/SheetTracker.jl`** — the
    existing PIII/PV/PVI RHS factories this layer wraps; their
    docstrings carry the FFW 2017 §2.1/§2.2 ground truth.
  - **`src/Problems.jl`** — `PadeTaylorProblem`, the struct this layer
    builds.
  - **Bead `padetaylor-9xf`** — the `path_network_solve` / `prob.order`
    footgun the forwarding methods must work around.
  - **Worklog 025** — where that footgun was caught.
  - Per-equation canonical forms: FW 2011 (PI), FW 2014 + FW 2015
    (PII), Reeger–Fornberg 2014 (PIV), FFW 2017 (PIII/PV/PVI) — all
    under `references/markdown/`.
