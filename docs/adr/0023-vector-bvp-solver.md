# ADR-0023 — First-order vector-system BVP solver (`VectorBVP.jl`)

**Status**: Accepted (2026-05-21) | **Bead**: `padetaylor-0ln.24` (VB1 —
`VectorBVP.jl` Chebyshev-collocation vector BVP solver) | **Plan**:
`docs/v0p2_plan.md` row VB1.

## Decision

`VectorBVP.vector_bvp_solve(f, z_a, z_b, B_a, B_b, g; …)` solves a
first-order **vector** ODE system on a complex segment

    y'(z) = f(z, y),    y : [z_a, z_b] ⊂ ℂ → ℂ^d,

subject to a **general linear two-point boundary condition**

    B_a · y(z_a) + B_b · y(z_b) = g,

`d` scalar equations with `B_a, B_b ∈ ℂ^{d×d}` and `g ∈ ℂ^d`. The method
is **Chebyshev spectral collocation + Newton**, the vector first-order
analogue of the existing scalar second-order `BVP.bvp_solve`.

Concretely the design fixes:

- **Discretisation.** `N+1` Chebyshev extrema nodes `t_j = cos(jπ/N)` on
  `[-1,1]` (DMSUITE convention, `t_0 = +1`, `t_N = -1`), affine-mapped to
  the segment by `z(t) = (z_a+z_b)/2 + s·t`, `s = (z_b-z_a)/2`. In `t`
  the system becomes `dY/dt = s·f(z(t), Y)`. The first-derivative matrix
  `D1` is built by the Trefethen/Weideman–Reddy formula — the construction
  is **copied verbatim** from `BVP._chebyshev_D1` (the private symbol is
  not reached into; the ~20-line builder is duplicated, per CLAUDE.md
  Rule "do NOT reach into a private symbol of another module").

- **Unknown layout.** The state is the `d`-vector `Y_j` at each of the
  `N+1` nodes, solved as a single stacked vector of length `(N+1)·d`,
  node-major (`[Y_0; Y_1; …; Y_N]`). The collocation operator on the
  stack is `D1 ⊗ I_d`.

- **τ-method row replacement.** Collocation gives `(N+1)·d` residual
  equations; the BC adds `d` more — the system is over-determined by `d`.
  The `d` collocation rows at node `t_0` (the `z_a` endpoint) are
  **replaced** by the `d` BC rows `B_a·Y_0 + B_b·Y_N - g`. This is the
  standard τ-method boundary-row substitution (Trefethen SMIM ch. 13:
  "replace rows of the differentiation matrix by the boundary condition").
  Replacing the `t_0` rows (rather than `t_N`, or splitting) is the
  documented choice; it is verified by the linear oracle reaching machine
  precision (test `VB.1.x`). Had that oracle floored above `~1e-12` the
  row choice would have been wrong and re-investigated per CLAUDE.md
  Rule 2 — it did not.

- **Newton with step-norm convergence.** Solve `J·ΔY = R`, update
  `Y ← Y - ΔY`, stop at `‖ΔY‖_∞ ≤ tol`. Step-norm (not residual-norm) is
  the convergence test for the same reason as scalar `BVP.bvp_solve`: the
  discrete residual floors at `cond(D1⊗I)·eps`, unachievable as a Newton
  target. The default `tol = eps(real(T))^(3/4)`.

- **Autodiff Jacobian via `Taylor1`, with optional override.** The `d×d`
  node Jacobian `∂f/∂y(z_j, Y_j)` is obtained without a new dependency:
  the `d`-vector `Y_j` is promoted to `Vector{Taylor1{CT}}` (constant
  series); for column `i` component `i` is seeded with a unit linear
  series, `f` is evaluated, and coefficient `[1]` of each output component
  is read off — `d` RHS evaluations give the full `d×d` matrix. This is
  exact (not finite-difference) and rides on the element-type-generic RHS
  convention of ADR-0020. A caller may pass `jacobian = Jf`, a function
  `(z,y) -> d×d matrix`, to bypass autodiff (cheaper for hand-derived
  Jacobians, and an escape hatch for an RHS that is not `Taylor1`-clean).

- **Barycentric vector callable.** Dense output is Berrut–Trefethen 2004
  Chebyshev-2 barycentric interpolation applied component-wise; the
  solution object is callable, `sol(z) -> Vector{CT}` of the `d`
  components. Out-of-segment `z` throws `DomainError` (mirrors `BVP.jl`).

- **Generic in `T`.** Float64 routes the Newton linear solve through
  LAPACK `\`; BigFloat routes through `GenericLinearAlgebra` (the same
  idiom as `BVP.jl`). `Complex{T}` endpoints propagate through `s`.

## Context

The motivating target is the `P_I⁽²⁾` tritronquée solution (bead
`padetaylor-r2m`; KKG 2015 §7). The tritronquée is a **separatrix**:
forward IVP integration diverges because the wrong-direction mode grows
exponentially, so a global (collocation / BVP) method is **mandatory** —
shooting cannot pin a separatrix. The v0.2 vector pipeline therefore
needs a BVP solver for first-order vector systems.

The existing scalar `BVP.bvp_solve` cannot be reused: it is hard-wired to
a **second-order scalar** ODE `u'' = f(z,u[,u'])` with **Dirichlet**
boundary conditions `u(z_a)=u_a, u(z_b)=u_b`. The `P_I⁽²⁾` companion form
is a first-order **4-vector** system, and the tritronquée boundary data
is naturally a general linear condition (asymptotic-amplitude pinning at
both ends), not pointwise Dirichlet on a scalar. A new module is the
honest answer.

## Alternatives considered, rejected

- **Reuse the scalar `BVP.jl` second-order path.** Rejected: the state,
  the `D2` operator, the Dirichlet BC enforcement, and the diagonal
  analytic Jacobian are all scalar-second-order specific. Lifting it in
  place would either break the v0.1 scalar tests or balloon the file past
  the 200-LOC cap. Additive — a new module — is the architecture rule.

- **Fourth-order scalar BVP.** The `P_I⁽²⁾` member is a fourth-order
  scalar ODE; one could build a scalar 4th-order Chebyshev BVP solver.
  Rejected: the v0.2 pipeline already represents `P_I⁽²⁾` as a first-order
  4-vector companion system (ADR-0022, `PainleveHierarchy.jl`); a
  first-order vector solver composes with that representation directly and
  generalises to the Noumi–Yamada `A_n⁽¹⁾` vector systems for free.

- **Shooting (IVP + Newton on the missing initial data).** Rejected for
  the headline reason above: the tritronquée is a separatrix, so any
  forward shoot diverges before reaching `z_b`. Collocation is global and
  does not integrate forward, so it is immune to the separatrix
  instability. Shooting is the wrong tool here by construction.

- **Finite-difference Jacobian.** Rejected: the `Taylor1` autodiff
  Jacobian is exact to machine precision and the RHS is already
  element-type-generic (ADR-0020), so autodiff is free and strictly more
  accurate than a finite-difference stencil.

## Consequences

- A new module `src/VectorBVP.jl` (literate top docstring; targeted at the
  Rule 6 ≤200-effective-LOC cap) and its test file
  `test/vector_bvp_test.jl` (prefix `VB.x.y`), with an inline
  mutation-proof block. Additive: no existing `src/` module is modified;
  only `test/runtests.jl` gains an `include`.
- The public API `vector_bvp_solve(f, z_a, z_b, B_a, B_b, g; N, tol,
  maxiter, initial_guess, jacobian) -> VectorBVPSolution` is the stable
  surface that VB2/VB3/V8b will depend on; `d` is inferred from
  `length(g)`.
- Wiring `VectorBVP` into the `PadeTaylor` umbrella module was deferred to
  VB3 — the VB1 test file `include`s `src/VectorBVP.jl` directly.  **Done
  in VB3** (bead `padetaylor-0ln.26`): `VectorBVP` is `include`d in
  `src/PadeTaylor.jl` immediately before `PainleveHierarchy` (so the latter
  can add a `PainleveHierarchyProblem` method of `vector_bvp_solve`), and
  `vector_bvp_solve` / `VectorBVPSolution` are re-exported.  VB3 also adds
  `PainleveHierarchy.painleve_hierarchy_jacobian(:I, m; t)` — the exact
  analytic companion-form Jacobian `∂f/∂y` of the `P_I^(m)` RHS — and a
  `vector_bvp_solve(php::PainleveHierarchyProblem, B_a, B_b, g; …)`
  convenience that boundary-value-solves a hierarchy member, defaulting
  `jacobian` to that analytic helper.  The helper is cross-validated
  against the `Taylor1` autodiff Jacobian and a finite-difference
  Jacobian (test `VW.2`).
- Fail-loud (Rule 1): `N < 4`, `maxiter < 1`, `tol < 0`, `B_a`/`B_b`/`g`/
  RHS dimension mismatches, Newton non-convergence, a singular Jacobian,
  and out-of-segment evaluation all throw with a `suggestion`/`detail`.

## References

- `src/BVP.jl` — the scalar second-order Chebyshev-Newton BVP solver this
  module mirrors structurally; `BVP._chebyshev_D1` is the verbatim source
  of the `D1` builder; the generic-`T` linear-solve and barycentric
  conventions are inherited from it.
- `references/bvp_recipe.md` — §1 Chebyshev nodes + affine map, §2 the
  `D1` matrix, §5 barycentric evaluation (the canonical recipe).
- `references/markdown/TrefethenSMIM_2000_book/…` — SMIM ch. 6 (`D1`),
  ch. 13 (Newton on nonlinear BVP, τ-method boundary-row replacement).
- `references/markdown/BerrutTrefethen2004_barycentric_SIAMReview/…` —
  the Chebyshev-2 barycentric formula and weights.
- ADR-0020 — the element-type-generic vector RHS convention the autodiff
  Jacobian relies on.
- ADR-0022 — the `P_I⁽²⁾` first-order vector companion-form representation.
- ADR-0001 (four-layer architecture), ADR-0002 (BigFloat linear algebra
  via GenericLinearAlgebra), ADR-0004 (Tier-3 BVP plan).
- CLAUDE.md Rules 1, 2, 4, 5, 6, 9, 10.
