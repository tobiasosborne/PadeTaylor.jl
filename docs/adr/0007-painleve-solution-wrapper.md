# ADR-0007 — `PainleveSolution` self-describing solve-output wrapper

**Status**: Accepted (2026-05-14) | **Bead**: `padetaylor-p1b`
**Context**: Design discussion 2026-05-14 — "zoom out and build some
nice Painlevé infrastructure; develop this into a useful package."
The first foundational piece: a self-describing solution object that
the rest of the Painlevé layer (named transcendents, a unified
`solve`, plotting) can hang off.

## Context

`PadeTaylorProblem` construction for the six Painlevé equations was
unified by ADR-0006's `PainleveProblem` builder.  But the *output*
side is still fragmented.  The solve layer produces **six unrelated
solution types**, each with its own shape and conventions:

| solver                        | output type            | shape                         |
|-------------------------------|------------------------|-------------------------------|
| `solve_pade`                  | `PadeTaylorSolution`   | callable dense trajectory     |
| `path_network_solve`          | `PathNetworkSolution`  | visited tree + grid scatter   |
| `bvp_solve`                   | `BVPSolution`          | callable Chebyshev interpolant|
| `dispatch_solve`              | `DispatcherSolution`   | 1-D IVP↔BVP chain             |
| `lattice_dispatch_solve`      | `LatticeSolution`      | 2-D lattice                   |
| `edge_gated_pole_field_solve` | `EdgeGatedSolution`    | masked pole-field grid        |

None of them is self-describing — a `PadeTaylorSolution` returned from
`solve_pade(::PainleveProblem)` carries no record that it *is* a
Painlevé-I solution, with which parameters, on which solution branch.
None pretty-prints.  `poles(...)` exists for exactly one of them
(`extract_poles(::PathNetworkSolution)`).  And ADR-0006's
`:transformed` forwarding methods (PIII / PV / PVI) simply **throw** —
a user who builds a transformed-frame problem cannot solve it through
the wrapper at all.

A "useful Painlevé package" needs a single solution abstraction that
(a) records its own provenance, (b) reads the same way regardless of
which solver produced it, (c) maps transformed-frame results back to
the natural `z`-frame, and (d) gives downstream layers — named
special-solution constructors, plotting — one type to target.

## Decision

Add `PainleveSolution`, a thin **self-describing wrapper** around
whichever raw solve-output was produced — the output-side mirror of
`PainleveProblem`.

### The wrapper struct

```julia
struct PainleveSolution{S, TF, FF}
    equation   :: Symbol                 # :I … :VI
    params     :: NamedTuple             # (), (α=…,), (α=…,β=…,γ=…,δ=…)
    name       :: Union{Symbol, Nothing} # :tritronquee, … or nothing (generic IC)
    frame      :: Symbol                 # :direct | :transformed
    to_frame   :: TF                     # (z,u,up) → solve-frame  (identity if :direct)
    from_frame :: FF                     # solve-frame → (z,u,up)  (identity if :direct)
    raw        :: S                      # the underlying solve output
end
```

`to_frame` / `from_frame` are carried verbatim from the
`PainleveProblem` the solution came from — the same `(z,u,up) ↔
(ζ,w,wp)` triple maps.  `name` is `nothing` for a problem built from
explicit ICs; it is reserved for the deferred named-transcendent
constructors (ADR-0006 "Deferred"), which will set it to
`:tritronquee`, `:hastings_mcleod`, etc.  Storing it now means that
later thread is purely additive — no struct change.

### What the wrapper provides

  1. **`show`** — a multi-line provenance summary: the equation's
     common name, its parameters, the named solution (if any), the
     frame, the `z`-frame domain, and a structural description of the
     wrapped raw output (segment count / visited-node count).  `show`
     stays cheap and total — it never computes poles (extraction can
     be expensive and can fail loud).

  2. **Callable `sol(z)`** — dense evaluation in the natural
     `z`-frame.  For a `:direct` problem whose `raw` is callable
     (`PadeTaylorSolution`) it forwards transparently.  For a
     `:transformed` problem it maps `z → ζ` via `to_frame`, evaluates
     `raw(ζ)`, and maps the result back via `from_frame` — see
     "The `:transformed` flip" below.  For a `raw` that is **not** a
     dense trajectory (`PathNetworkSolution` — a grid scatter) it
     **throws** with a suggestion (CLAUDE.md Rule 1): a grid solve has
     no dense interpolant, and silently returning a nearest-grid value
     would be a truthful-looking lie.

  3. **`poles(sol)`** — pole locations in the `z`-frame.  Forwards to
     `PoleField.extract_poles` for a `PathNetworkSolution` raw, and
     (via bead `padetaylor-26r`) for a `PadeTaylorSolution` raw.
     For a `:transformed` problem the extracted `ζ`-frame poles are
     mapped back through `from_frame`.

  4. **`grid_values(sol)`** — `(z, u, u')` point values in the
     `z`-frame.  For a `PathNetworkSolution` raw this is the Stage-2
     grid; for a `PadeTaylorSolution` raw it is the segment
     breakpoints.  For a `:transformed` problem each `(ζ, w, w')`
     triple is mapped back through `from_frame` — this is what makes a
     `:transformed` path-network result actually consumable (its
     `raw.grid_*` are `ζ`-frame).

  5. **Accessors** — `equation(sol)`, `parameters(sol)`,
     `solutionname(sol)` — read the provenance without reaching into
     fields.

### The `:transformed` flip

ADR-0006 refinement #2 decided that `:transformed` forwarding methods
*throw*, because a `ζ`-frame `PathNetworkSolution` "cannot be
faithfully round-tripped — its per-node Padé store and step length
`h` are intrinsically `ζ`-frame and have no `z`-frame image."

That reasoning is correct **about the Padé store** and stays correct.
But it conflated two things: the *internal* `ζ`-frame machinery (Padé
coefficients, `h`), which genuinely has no `z`-image, and the *point
values* (`grid_u`, `visited_u`, the callable's `(w, wp)` output),
which map back cleanly via `from_frame`.  A solution object that
exposes point values in the `z`-frame while keeping its Padé store
`ζ`-frame is honest, not half-mapped — the user asked for `u(z)`, gets
`u(z)`; the `ζ`-frame internals are an implementation detail they
never see.

`PainleveSolution` is exactly the abstraction that makes this clean:
it carries `to_frame` / `from_frame`, so the callable, `poles`, and
`grid_values` do the mapping at the boundary.  Therefore:

  - `path_network_solve(pp::PainleveProblem, grid)` **no longer throws**
    for `:transformed` problems.  It maps the caller's `z`-frame `grid`
    into the `ζ`-frame, solves there, and returns a `PainleveSolution`
    whose `poles` / `grid_values` present `z`-frame values.

  - `solve_pade(pp::PainleveProblem)` is a **partial** flip, and the
    reason is honest: `solve_pade` does fixed-step *real-axis*
    stepping, and a `:transformed` problem's `ζ`-domain is in general
    *complex* (`ζ = log z`), so `state.z < z_end` is not even defined.
    `solve_pade` therefore serves a `:transformed` problem only in the
    special case where the `ζ`-domain is real-typed (the IC and span
    were given as real numbers); for the common complex-`ζ` case it
    **throws** — but now with a message that points the caller at
    `path_network_solve`, not the old "cannot be round-tripped".

  - The `ζ`-frame Padé store remains reachable as `sol.raw` for a
    caller who genuinely wants the transformed-frame internals; the
    docstring says so plainly.  `sol.raw.grid_z` for a `:transformed`
    path-network result is `ζ`-frame — `grid_values(sol)` is the
    `z`-frame view.

This supersedes ADR-0006 refinement #2.  The deferred bead
`padetaylor-soi` ("`:transformed` forwarding with point-value
remapping") is **subsumed** by this ADR and closed.

### Packaging — one module, two files

`PainleveSolution` lives in `src/PainleveSolution.jl` but in the
**same `Painleve` module** as `PainleveProblem`: `src/Painleve.jl`'s
module body does `include("PainleveSolution.jl")` before the
constructor / forwarding code that needs the type.  `PainleveProblem`
and `PainleveSolution` are one conceptual namespace and belong
together; CLAUDE.md Rule 6's ≤200 LOC budget is **per file**, and the
two-file split keeps each well under it (the alternative — a separate
`PainleveSolutions` module — would force an awkward
problem→solution module dependency for no conceptual gain).

## Scope

**v1 (bead `padetaylor-p1b`):**
  - `src/PainleveSolution.jl` — the struct, `show`, callable, `poles`,
    `grid_values`, accessors.
  - `src/Painleve.jl` — `include` the new file; rewrite the
    `solve_pade` / `path_network_solve` forwarding methods to wrap
    their result in a `PainleveSolution`.  `path_network_solve` flips
    fully for `:transformed`; `solve_pade` flips for the real-`ζ`
    case and throws a `path_network_solve`-pointing message otherwise.
  - Re-export `PainleveSolution`, `poles`, `grid_values`, `equation`,
    `parameters`, `solutionname` from the umbrella module.
  - Tests: provenance fields populated correctly per equation; the
    frame-mapped callable matches a direct
    `from_frame(raw(to_frame(z)…))` computation; fail-loud on a
    grid-type callable; `poles` round-trips on a `:direct`
    path-network solution.  Mutation-prove (perturb the frame map →
    RED; drop the grid-type guard → RED).

**v1, separate beads (depend on `padetaylor-p1b`):**
  - `padetaylor-26r` — extend `PoleField.extract_poles` to accept a
    `PadeTaylorSolution` (single-trajectory pole extraction), so
    `poles(sol)` works on `:direct` IVP results, not only
    path-networks.
  - `padetaylor-ylr` — a Makie plot recipe package extension
    (`ext/PadeTaylorMakieExt.jl`) targeting `PainleveSolution`.

**Deferred (not v1):**
  - `PainleveSolution` constructors / forwarding for the other four
    output types (`BVPSolution`, `DispatcherSolution`,
    `LatticeSolution`, `EdgeGatedSolution`) — the struct's `S`
    parameter is already generic over them; what is missing is
    `PainleveProblem` forwarding methods for those solvers, which is a
    separate piece of work.
  - Named special-solution constructors that populate `name`
    (`tritronquee(:I)`, `hastings_mcleod()`, …) — ADR-0006 "Deferred".
  - A unified `solve(pp; over=…)` that picks IVP / path-network / BVP
    from the requested domain.

## Fail-loud behaviour (CLAUDE.md Rule 1)

  - Calling `sol(z)` on a `PainleveSolution` whose `raw` is a grid
    type (`PathNetworkSolution`) throws an `ArgumentError`: a
    path-network is a grid of independent evaluations, not a dense
    interpolant; suggest reading `sol.raw.grid_*` directly or solving
    with `solve_pade` for a callable trajectory.
  - `poles(sol)` on a `raw` type with no pole-extraction route throws
    rather than returning an empty vector — an empty result must mean
    "no poles found", never "extraction not wired".
  - Evaluating a `:transformed` callable at a `z` the transform
    cannot reach (a fixed branch point) propagates the underlying
    transform's `DomainError` unchanged.

## Rejected alternatives

  1. **Make every raw solve-output self-describing in place** (add
     `equation` / `params` fields to `PadeTaylorSolution` etc.) —
     rejected: those types are equation-agnostic by design (ADR-0001's
     solver substrate takes an arbitrary `f(z,u,u')`); bolting
     Painlevé identity onto them inverts the dependency.  The wrapper
     keeps the substrate clean.

  2. **A separate `PainleveSolutions` module** — rejected: it would
     force `Painleve` (the problem builder) to depend on it just so
     the forwarding methods can construct the type, for no conceptual
     gain.  One module, two files is cleaner.

  3. **Keep `:transformed` forwarding throwing** (leave ADR-0006
     refinement #2 standing, do point-value remapping in a later
     bead) — rejected: the wrapper *is* the remapping vehicle;
     building it and then not using it for the transformed equations
     would ship the package with half its equations still
     second-class.  The honest move is to flip the decision now and
     document that `sol.raw` is the `ζ`-frame escape hatch.

  4. **Subtype `PainleveSolution <: PadeTaylorSolution`** — rejected
     for the same reason ADR-0006 rejected `PainleveProblem <:
     PadeTaylorProblem`: the wrapper must wrap *six* different raw
     types, not specialise one.  A generic `{S}` wrapper with
     forwarding methods is the right shape.

## Consequences

  - **Positive**: one self-describing solution type for the whole
    Painlevé layer; `solve_pade(pp)` / `path_network_solve(pp)` return
    something that prints its own provenance and knows how to find its
    poles.  Downstream threads (named transcendents, plotting, a
    unified `solve`) all target one type.
  - **Positive**: PIII / PV / PVI become first-class — the
    `:transformed` equations can finally be solved end-to-end through
    the wrapper, with `z`-frame point values out.
  - **Positive**: thin and additive — no change to the six raw
    solution types or the solver substrate; the existing test suite is
    unaffected except for the new `PainleveSolution` tests and the
    rewritten forwarding-method tests.
  - **Negative / honest**: the wrapper does not unify the *internal*
    representations — a `PainleveSolution` over a `PadeTaylorSolution`
    and one over a `PathNetworkSolution` still behave differently
    (callable vs. not).  The wrapper unifies *provenance and access*,
    not *representation*.  This is the right boundary: forcing a
    common internal representation on a dense trajectory and a grid
    scatter would lose information either way.
  - **Negative / honest**: for `:transformed` problems the Padé store
    in `sol.raw` is `ζ`-frame and stays that way.  A caller who
    reaches past the wrapper into `sol.raw.pade` / `sol.raw.h` /
    `sol.raw.visited_*` gets `ζ`-frame numbers.  The docstring flags
    this; the supported surface (`sol(z)`, `poles(sol)`) is fully
    `z`-frame.

## References

  - **ADR-0006** — the `PainleveProblem` builder this mirrors;
    refinement #2 (the `:transformed`-throws decision) is superseded
    here.
  - **ADR-0001** — the four-layer architecture; `PainleveSolution`
    sits beside `PainleveProblem` as an output-side convenience above
    the solver substrate.
  - **ADR-0004** — the path-network's per-node Padé store that
    `extract_poles` (and hence `poles`) reads.
  - **CLAUDE.md Rule 1** — fail loud; the grid-type callable guard and
    the `poles`-not-wired guard both derive from it.
  - **CLAUDE.md Rule 6** — ≤200 LOC per file; the motivation for the
    two-file / one-module split.
  - **`src/Problems.jl`**, **`src/PathNetwork.jl`** — the
    `PadeTaylorSolution` / `PathNetworkSolution` raw types wrapped.
  - **`src/CoordTransforms.jl`**, **`src/SheetTracker.jl`** — the
    `z ↔ ζ` transform pairs carried as `to_frame` / `from_frame`.
  - **Bead `padetaylor-soi`** — subsumed and closed by this ADR.
