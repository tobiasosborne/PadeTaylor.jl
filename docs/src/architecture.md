# Architecture

PadeTaylor.jl solves Painlevé-type ODEs in the complex plane using a
Taylor–Padé integration strategy derived from Fornberg & Weideman 2011
(FW 2011) and the robust Padé algorithm of Gonnet, Güttel & Trefethen
2013 (GGT 2013).  This chapter explains how the package is layered, why
each architectural choice was made, and what the separation of concerns
buys you as either a user or a contributor.

Architecture-decision records that this chapter synthesises live under
`docs/adr/` in the repository:

  - **ADR-0001** — `docs/adr/0001-four-layer-architecture.md`
  - **ADR-0002** — `docs/adr/0002-bigfloat-svd-via-genericlinalg.md`
  - **ADR-0003** — `docs/adr/0003-extensions-pattern.md`
  - **ADR-0004** — `docs/adr/0004-path-network-architecture.md`

## Four-layer architecture (ADR-0001)

The core algorithm decomposes into four algorithmically independent
layers, each a separate source module:

```
1. LinAlg / Coefficients  ─── compute c₀, …, cₙ of the local Taylor
                               expansion of y(z + h) at state (z, y).
2. RobustPade             ─── convert the truncated series to a diagonal
                               (n/2, n/2) Padé approximant via GGT 2013
                               Algorithm 2 + Chebfun reweighting.
3. StepControl            ─── choose step length h from coefficient decay
                               (Jorba–Zou 2005) or denominator-root
                               distance (FW 2011 §3.1 path heuristic).
4. PadeStepper            ─── orchestrate one step: coefficients → Padé →
                               step length → evaluate at the new point.
```

This decomposition mirrors FW 2011's own exposition: §2.1 (Taylor
coefficients), §2.2 (Padé conversion), §3.1 (path-direction selection),
§5.1 (order choice).  **The layers are the natural fault lines of the
algorithm.**  A bug in the coefficient layer is visible before any Padé
conversion — checkable against analytic derivatives.  A bug in the Padé
layer surfaces on known-correct input coefficients — checkable against
`padeapprox.m`.  A bug in step control surfaces on known-correct
coefficients and Padé output — checkable against Jorba–Zou's step
formula directly.  A bug in the orchestrator surfaces as wrong
end-of-step state — checkable against FW Table 5.1.

The interface between adjacent layers is a typed Julia struct
(`PadeApproximant{T}`, `PadeStepperState{T}`), not a positional tuple.
Caller code never depends on the internal layout of a lower layer.

Above these four layers sits the **driver tier** (`Problems.solve_pade`,
and above that `PathNetwork.path_network_solve`).  The driver tier
composes the inner layers but does not duplicate their logic: every
re-entry point for new ODE families (`BVP.bvp_solve`,
`Dispatcher.dispatch_solve`, `LatticeDispatcher.lattice_dispatch_solve`)
calls the same `PadeStepper.pade_step_with_pade!` primitive and inherits
the same coefficient, Padé, and step-control modules unchanged.

## BigFloat-SVD via `GenericLinearAlgebra` (ADR-0002)

GGT 2013 Algorithm 2 classifies a singular value as "zero" when
`σᵢ < tol · ‖c‖₂`.  Near a block boundary in the Padé table the gap
between the smallest genuine singular value and the largest noise
singular value can be narrow — far narrower than the dynamic range
Julia's LAPACK-backed `LinearAlgebra.svd` guarantees for `Float64`.

The SVD dispatcher in `LinAlg.pade_svd` resolves this with a
precision-indexed dispatch table:

| Element type | Backend | Guarantee |
|---|---|---|
| `Float64` / `Float32` | `LinearAlgebra.svd` (LAPACK) | Matches `padeapprox.m` exactly; absolute-error Demmel–Kahan suffices at Float64. |
| `BigFloat` | `GenericLinearAlgebra.svd` | One-sided Jacobi (Demmel–Veselić): *relative* error `c · 2⁻ᵖ · σᵢ` per singular value — small-but-genuine SVs stay above the GGT threshold. |
| `Arblib.Arb` | Extension: convert to `Matrix{BigFloat}`, then as above | `Arblib.jl` has no SVD (verified by source inspection; see `RESEARCH.md §5.1`). |

For GGT matrices at the typical sizes (`n ≤ 60`) the runtime overhead
of Jacobi over Demmel–Kahan is negligible.
`GenericLinearAlgebra.jl` is chosen because it is the only
production-grade Julia library providing generic SVD over
`T <: AbstractFloat`, including `BigFloat`, with no FFI surface.

**One precision-loss caveat applies to the `Arb` path.**  Converting
`Arb(mid ± rad)` to `BigFloat` discards the ball radius.  This is
acceptable because the downstream consumer of the SVD is
`RobustPade.robust_pade`, which uses the results to build a
`PadeApproximant{T}` whose precision rests on the GGT 2013
normalisation `‖b‖₂ = 1`.  The reported precision of the Padé
approximant is set by the *coefficient* arithmetic (where Arb radii
are correctly tracked), not by the SVD step itself.  A fully rigorous
Arb-arithmetic SVD would require a verified-arithmetic Jacobi port;
that is explicitly deferred beyond v1 scope.

## Extensions pattern (ADR-0003)

`using PadeTaylor` loads three hard dependencies — `TaylorSeries`,
`GenericLinearAlgebra`, `Polynomials` — all pure-Julia, no FFI, no
large precompile cost.  Two optional features live in Julia 1.9+
*package extensions* (weak dependencies):

```toml
[weakdeps]
Arblib      = "..."
CommonSolve = "..."

[extensions]
PadeTaylorArblibExt      = "Arblib"
PadeTaylorCommonSolveExt = "CommonSolve"
```

**`PadeTaylorArblibExt`** (`ext/PadeTaylorArblibExt.jl`) activates when
the user also loads `Arblib`.  It adds the `Matrix{Arb} → BigFloat`
shim in `LinAlg.pade_svd` described above, plus an
`Arb`-precision-aware `default_tol`.  A user who only needs `Float64`
integration pays no precompile cost for Flint's C-FFI surface.

**`PadeTaylorCommonSolveExt`** (`ext/PadeTaylorCommonSolveExt.jl`)
activates when the user loads `CommonSolve`.  It provides the
`PadeTaylorAlg <: CommonSolve.AbstractAlgorithm` marker and implements
the `init` / `step!` / `solve!` interface so PadeTaylor.jl participates
in the wider SciML solver ecosystem.  The extension contains no
algorithmic logic — it is a translation layer only, delegating every
call to `Problems.solve_pade`.

Both extensions are loaded automatically by Julia's extension mechanism;
the user does not import anything extra beyond `using PadeTaylor; using
Arblib` (or `using CommonSolve`).

## Path-network architecture (ADR-0004)

Naïve real-axis fixed-step integration cannot reach FW Table 5.1
accuracy: the Painlevé singularities cluster along the real axis and
step-size heuristics cannot reliably vault them.  FW 2011 §3.1 provides
the remedy — a **5-direction complex-plane path-network** that builds a
visited-node tree in the complex plane and then fills a fine output
grid by barycentric extrapolation from stored Padé approximants.

`PathNetwork.path_network_solve` implements this in two stages.

**Stage 1** (FW 2011 §3.1,
`FW2011_painleve_methodology_JCP230.md:155–166`) builds the
visited-node tree.  Starting from the initial condition, each requested
output point is reached by launching up to five candidate steps in the
angular wedge `[-π/4, -π/8, 0, π/8, π/4]` (the FW 2011 default;
overridable via the `wedge_angles` keyword) and selecting the candidate
that minimises `|u|` (the `:min_u` heuristic, or optionally
`:steepest_descent` via `θ = arg(-u/u')`).  Each successful intermediate
node is stored in the visited set together with its full
`PadeApproximant` — this is the **canonical-Padé invariant**: the
stored Padé evaluated at `t = 0` must recover the stored solution
value exactly.

**Stage 2** (`FW2011_painleve_methodology_JCP230.md:166, 397`) fills
the fine output grid without any new Taylor jets.  For each fine-grid
point `z_f`, the algorithm finds the nearest Stage-1 visited node
`z_v` with coverage `|z_f - z_v| ≤ h`, retrieves its stored Padé, and
evaluates at the local coordinate `t = (z_f - z_v) / h`.  If no
visited node is within coverage, the output is set to `NaN` with a
diagnostic; silent extrapolation is never performed.

The `path_network_solve` signature is:

```julia
path_network_solve(
    prob             :: PadeTaylorProblem,
    grid             :: AbstractVector{<:Complex};
    h                :: Real    = 0.5,
    order            :: Integer = 30,
    wedge_angles     :: AbstractVector{<:Real} = [-π/4, -π/8, 0.0, π/8, π/4],
    step_selection   :: Symbol  = :min_u,
    step_size_policy :: Symbol  = :fixed,
    max_steps_per_target :: Integer = 1000,
    rtol             :: Real    = 1e-10,
) -> PathNetworkSolution
```

The `PathNetworkSolution` struct carries both the Stage-1 visited tree
(`visited_z`, `visited_u`, `visited_pade`) and the Stage-2 interpolated
grid (`grid_z`, `grid_u`, `grid_up`), so callers can inspect the
intermediate path structure for diagnostics.

### Composition layers above the path-network

Three further driver layers share the same inner primitives:

- **`BVP.bvp_solve`** — Chebyshev–Newton spectral BVP solver
  (FW 2011 §3.2, Trefethen *SMIM* ch. 6 & 13), sharing the
  `Coefficients` and `RobustPade` layers, replacing `StepControl` and
  `PadeStepper` with a spectral collocation Newton loop.

- **`Dispatcher.dispatch_solve`** — 1D IVP↔BVP chain composition:
  segments detected as smooth by `EdgeDetector` are solved by BVP
  fill; others by IVP stepping.

- **`LatticeDispatcher.lattice_dispatch_solve`** — 2D-grid
  composition: per-row BVP fill on smooth runs flanked by IVP cells,
  per FW 2011 line 190.

Coordinate-transform helpers for PIII and PV (`CoordTransforms`) and
the Riemann-sheet winding primitives for PVI (`SheetTracker`) wrap the
path-network externally — no changes to the inner layers are required.

## Where to read next

- **`LinAlg`** (`src/LinAlg.jl`) — SVD dispatcher.
- **`RobustPade`** (`src/RobustPade.jl`) — GGT 2013 Algorithm 2 +
  Chebfun reweighting; `robust_pade`, `PadeApproximant`.
- **`Coefficients`** (`src/Coefficients.jl`) — Taylor jet generation;
  `taylor_coefficients_1st`, `taylor_coefficients_2nd`.
- **`StepControl`** (`src/StepControl.jl`) — `step_jorba_zou`,
  `step_pade_root`.
- **`PadeStepper`** (`src/PadeStepper.jl`) — `pade_step!`,
  `pade_step_with_pade!`.
- **`Problems`** (`src/Problems.jl`) — `PadeTaylorProblem`,
  `solve_pade`, `PadeTaylorSolution`.
- **`PathNetwork`** (`src/PathNetwork.jl`) — `path_network_solve`,
  `PathNetworkSolution`.
- **`BVP`** (`src/BVP.jl`) — `bvp_solve`.
- **`Dispatcher` / `LatticeDispatcher`** (`src/Dispatcher.jl`,
  `src/LatticeDispatcher.jl`) — chain and 2D-grid composition.
- **`CoordTransforms`** (`src/CoordTransforms.jl`) — PIII / PV
  exponential-coordinate helpers.
- **`SheetTracker`** (`src/SheetTracker.jl`) — PVI Riemann-sheet
  winding primitives.

Full per-module docstrings appear in the [API](api.md) reference.
