# ADR-0024 — Harmonic-extension sector fill: triple-method majority vote (`Laplace2D.jl`)

**Status**: Accepted (2026-05-21) | **Bead**: `padetaylor-0ln.35` (in-house
2D-Chebyshev Laplace solver + this ADR) | **Plan**: `docs/v0p2_plan.md`,
the `P_I⁽²⁾` tritronquée surface figure (KKG 2015 Figs 7.4/7.5).

## Decision

The pole-free-sector fill of the `P_I⁽²⁾` tritronquée surface figure is
computed **three algorithmically independent ways** and **majority-voted
per Cartesian grid point**:

  1. **Ray-fan `vector_bvp_solve`** — a fan of complex line segments
     emanating through the sector, each a first-order vector Chebyshev
     collocation BVP (ADR-0023, `src/VectorBVP.jl`). The direct vector
     solution along each ray.

  2. **In-house 2D-Chebyshev spectral Laplace solve** — THIS bead's
     `src/Laplace2D.jl`. `Re u` and `Im u` of the tritronquée are
     **harmonic** (the tritronquée `u` is analytic on its pole-free
     sector), so each satisfies `∇²φ = 0` on the sector interior and is
     fully determined by its boundary values. A single tensor-product
     Chebyshev Laplacian solve on a rectangle recovers the interior.

  3. **Gridap.jl FEM Laplace solve** — a Gridap.jl finite-element
     solution of the same Dirichlet Laplace problem, an algorithmically
     disjoint third voter. It is a **`figures/`-project dependency + a
     plain figure-helper file** (`figures/_kkg_pi2_gridap_helper.jl`),
     **NOT a package weak-dep** — see the *Correction* section below for
     the root cause of that placement.

The three interior fields are reconciled by a **per-grid-point majority
vote / agreement map** (assembled in F3, bead `padetaylor-0ln.33`):
points where the three methods agree to tolerance are accepted; the
agreement map is itself a figure-quality diagnostic. This mirrors the
project's standing triple-oracle discipline (CLAUDE.md Rule 4; cf. the
`SharedPade` Padé routine cross-validated against Chebfun's
`padeapprox.m`, a closed-form oracle, and a finite-difference oracle).
One numerical method can hide a systematic bug behind a plausible-looking
surface; three independent methods that agree cannot all share it.

## Correction (2026-05-21) — Gridap is a `figures/` dependency, not a package weak-dep

The original plan (and commit `cfe8b86`) shipped voter (3) as a
`PadeTaylorGridapExt` **package weak-dep extension** with
`Gridap = "0.18"` in `Project.toml` `[weakdeps]` / `[extensions]` /
`[compat]` (the ADR-0003 extensions pattern). **This broke the package
build.** `julia --project=. -e 'using Pkg; Pkg.test()'` failed at
environment resolution with `Unsatisfiable requirements detected for
package ForwardDiff`: Gridap 0.18 pins `ForwardDiff 0.10.x`, which is
**unsatisfiable** against the core package's `SpecialFunctions 2.7.2`
stack (it needs a newer `ForwardDiff`). The package test environment
became unresolvable — the whole suite could not run.

The root cause is architectural, not a version typo. Gridap is a
~160-dependency FEM package. Making it a *package* weak-dep puts its
entire conflict surface into the core package's resolvable set even
though Gridap here is **figure-only infrastructure** — one of three
Laplace solvers cross-checked in a single headline figure, exercised
only inside F3. The ADR-0003 weak-dep pattern is right for *small*,
presentation-only deps (Makie, Arblib, CommonSolve, DelaunayTriangulation);
it is the wrong tool for a heavy FEM stack with a deep transitive
conflict surface.

**The correction (bead `padetaylor-0ln.36`, corrective phase):**

- Gridap is **removed** from the package `Project.toml` (`[weakdeps]`,
  `[extensions]`, `[compat]`, `[extras]`, `[targets]`); the
  `ext/PadeTaylorGridapExt.jl` extension and the
  `laplace2d_solve_gridap` stub in `src/Laplace2D.jl` are deleted; the
  package `Manifest.toml` is restored to its pre-`cfe8b86` state.
  `Pkg.test()` resolves cleanly again.
- Gridap is **re-homed into the separate `figures/` project**: added to
  `figures/Project.toml` `[deps]` with `Gridap = "0.20"` in `[compat]`.
  Gridap **0.20** (0.20.7 at time of writing) dropped the `ForwardDiff
  0.10` pin and resolves cleanly in the figures env alongside CairoMakie
  + dev PadeTaylor on Julia 1.12.
- The FEM solver becomes a **plain figure-helper file**,
  `figures/_kkg_pi2_gridap_helper.jl` (the `figures/_kkg_pi2_helpers.jl`
  pattern — figure scripts `include` it), exposing
  `laplace2d_solve_gridap(...)`. The logic is ported verbatim from the
  deleted extension; it still returns a `PadeTaylor.Laplace2DSolution`
  so F3 treats the FEM and spectral voters uniformly.
- The Gridap verification is re-homed to
  `figures/test_kkg_pi2_gridap.jl`, a standalone script run under
  `--project=figures` (closed-form harmonic functions + agreement with
  `laplace2d_solve` to FEM tolerance; ~5e-5 FEM-vs-spectral agreement,
  O(h²) convergence confirmed). The main `test/` suite stays Gridap-free.

**The triple-method majority vote is unchanged.** Voters (1) ray-fan
BVP, (2) in-house 2D-Chebyshev, (3) Gridap FEM all still feed the
per-grid-point agreement map — that vote simply executes in F3 *within
the `figures/` environment*, where Gridap is available, rather than via
a package extension.

## Context

The headline figure (KKG 2015 §7, Figs 7.4/7.5) renders `Re u` / `Im u`
of the `P_I⁽²⁾` tritronquée over the complex plane. The tritronquée is a
separatrix solution: it is **pole-free on an angular sector** and
pole-rich elsewhere. On the pole-rich complement the vector path-network
(`VectorPathNetwork`, Stage-2 fill `VectorPathNetworkStage2`) does the
job. On the pole-free sector the solution is analytic, hence `Re u` and
`Im u` are **harmonic** — and a harmonic function on a domain is the
unique solution of the Dirichlet Laplace problem with the function's own
boundary values. KKG 2015 §7 themselves fill the sector this way,
citing Trefethen's spectral-Laplace program.

A pole-free *angular sector* `{ r₀ ≤ |x| ≤ r₁, θ₀ ≤ arg x ≤ θ₁ }` is
**conformally a rectangle** under the map `w = log x`: the radial
coordinate becomes `Re w ∈ [log r₀, log r₁]` and the angular coordinate
becomes `Im w ∈ [θ₀, θ₁]`. The Laplacian is conformally invariant —
`φ` harmonic in `x` ⟺ `φ ∘ exp` harmonic in `w` — so the sector
Laplace problem is exactly a *rectangle* Laplace problem in `w`. That
makes a tensor-product Chebyshev solver on a rectangle the natural
in-house method.

## The in-house 2D-Chebyshev solver (`src/Laplace2D.jl`)

`laplace2d_solve` solves the Dirichlet problem `∇²φ = 0` on a rectangle
`[a,b] × [c,d]`:

- **Discretisation.** Chebyshev extrema grids on each axis: `Nx+1` nodes
  in `x`, `Ny+1` in `y` (DMSUITE convention, `t₀ = +1`, `t_N = -1`),
  affine-mapped to `[a,b]` and `[c,d]`. The first-derivative matrix is
  built by the Trefethen/Weideman–Reddy formula — the `_chebyshev_D1`
  builder is **copied verbatim** from `BVP.jl` (its private symbol is
  not reached into; the ~20-line builder is duplicated, the same
  discipline ADR-0023 records for `VectorBVP.jl`). `D₂ = D₁·D₁`, scaled
  by the affine factor `(2/(b-a))²` resp. `(2/(d-c))²` (chain rule on the
  `t ↔ x` map; squared because it is a *second* derivative).

- **Tensor-product Laplacian.** With nodes in lexicographic order the
  discrete Laplacian on the *interior* nodes is the Kronecker sum
  `L = kron(D²x[int,int], I_y) + kron(I_x, D²y[int,int])` (Trefethen
  SMIM Program 16 / Program 23, `references/markdown/TrefethenSMIM_2000_book/
  TrefethenSMIM_2000_book.md:753–760`).

- **Dirichlet BC into the RHS.** `D²` is split into interior×interior and
  interior×boundary blocks; the boundary blocks act on the known edge
  data and move to the right-hand side — the 2D analogue of `BVP.jl`'s
  `bc_col` boundary-contribution column.

- **One linear solve.** `φ_int = L \ rhs`. The problem is **linear** —
  no Newton iteration (unlike `BVP.bvp_solve`). Generic in `T`: Float64
  routes `\` through LAPACK, BigFloat through `GenericLinearAlgebra` (the
  `BVP.jl` idiom).

- **Dense output.** The solution object carries `φ` on the Chebyshev
  tensor grid plus the grids; it is callable, `sol(x,y) -> φ`, via 2D
  barycentric interpolation (tensor product of the Berrut–Trefethen 2004
  Chebyshev-2 1D barycentric formula).

`Laplace2D` is **domain-agnostic**: it solves Laplace on a plain
rectangle and knows nothing of the `w = log x` sector map or the
`P_I⁽²⁾` boundary data. Composing the conformal map and supplying the
tritronquée boundary data (outer arc from `pI2_tritronquee_ic`'s
asymptotic series; sector-edge rays from `vector_bvp_solve`) is F3's job
(bead `padetaylor-0ln.33`).

## Alternatives considered, rejected

- **A single method (no vote).** Rejected: the figure is a headline
  deliverable; a lone numerical surface can carry a systematic error
  (a sign slip in the conformal metric, an off-by-one in a boundary
  index) that looks entirely plausible. The triple-oracle discipline
  (Rule 4) is exactly the guard against that.

- **Gridap.jl / Ferrite.jl FEM as the *sole* method.** Rejected as sole
  method: the sector maps conformally to a rectangle, and on a rectangle
  a tensor-product spectral solve converges geometrically with a handful
  of nodes — unstructured FEM is overkill and converges only
  algebraically. Gridap is, however, **kept as the independent third
  voter** precisely because it is algorithmically disjoint from the two
  spectral/collocation methods: a FEM bug and a spectral bug will not
  coincide.

- **Gridap as a package weak-dep extension** (the original `cfe8b86`
  plan). Rejected — see the *Correction* section: Gridap 0.18's
  `ForwardDiff 0.10.x` pin is unsatisfiable against the core package's
  `SpecialFunctions` stack, and a ~160-dependency FEM package does not
  belong in the core package's resolvable set for figure-only
  infrastructure. Gridap lives in the `figures/` project instead.

- **ApproxFun.jl / MethodOfLines.jl.** Rejected: ApproxFun's PDE support
  would duplicate the spectral method we already own via `_chebyshev_D1`;
  MethodOfLines is a finite-difference method-of-lines discretiser aimed
  at time-dependent PDEs, the wrong tool for a static elliptic solve.
  Neither adds an *independent* third voter beyond what Gridap provides.

- **Bilinear interpolation across the ray fan.** Rejected: bilinear (or
  any local polynomial) interpolation between rays is not a principled
  reconstruction of a *harmonic* function — it does not satisfy
  `∇²φ = 0` in the interior and would smear the figure. A harmonic
  function must be filled by a harmonic solve.

## Consequences

- A new module `src/Laplace2D.jl` (literate top docstring; targeted at
  the Rule 6 ≤200-effective-LOC cap) and its test file
  `test/laplace2d_test.jl` (prefix `L2D.x.y`), with an inline
  mutation-proof block. `Laplace2D` is `include`d in `src/PadeTaylor.jl`
  and `laplace2d_solve` / `Laplace2DSolution` are re-exported.

- Bead `padetaylor-0ln.36` adds the Gridap FEM voter (third voter) as a
  `figures/` figure-helper (`figures/_kkg_pi2_gridap_helper.jl`) with
  `Gridap` in `figures/Project.toml` `[deps]` (`[compat]` `0.20`) and a
  Law-2 justification line — NOT a package weak-dep (see the *Correction*
  section). Its verification is `figures/test_kkg_pi2_gridap.jl`, run
  under `--project=figures`.

- F3 (bead `padetaylor-0ln.33`) composes the three methods *within the
  `figures/` environment*: it applies the `w = log x` conformal map,
  supplies the tritronquée boundary data, runs all three solvers
  (`vector_bvp_solve`, `laplace2d_solve`, the
  `_kkg_pi2_gridap_helper.jl` `laplace2d_solve_gridap`), and assembles
  the per-grid-point majority vote and the agreement map.

- Fail-loud (Rule 1): `Nx`/`Ny` too small, and `bc_*` boundary-data
  dimension mismatches, all throw with a `suggestion`/`detail`.

## References

- `references/markdown/TrefethenSMIM_2000_book/TrefethenSMIM_2000_book.md`
  — SMIM ch. 6 (`D₁`/`D₂`), ch. 7 Program 16 + ch. 9 Program 23
  (`:753–760`: the tensor-product Chebyshev Laplacian
  `L = kron(I,D²) + kron(D²,I)` on interior nodes, Dirichlet BC handling).
- `src/BVP.jl` — `_chebyshev_D1` (the verbatim source of the `D₁`
  builder, `:525–547`); the `D₂ = D₁²` + affine-scale and the
  `bc_col` boundary-contribution idiom (`:279–297`); the generic-`T`
  LAPACK/`GenericLinearAlgebra` linear-solve convention.
- `references/markdown/BerrutTrefethen2004_barycentric_SIAMReview/…` —
  the Chebyshev-2 barycentric formula and weights (the 2D callable is
  its tensor product).
- ADR-0023 — `VectorBVP.jl`, voter (1); also the precedent for copying
  `_chebyshev_D1` verbatim rather than reaching into a private symbol.
- ADR-0003 — the weak-dep extension pattern; the *Correction* section
  records why Gridap voter (3) is **not** placed there (heavy FEM stack
  + an unsatisfiable `ForwardDiff` pin) and lives in `figures/` instead.
- `figures/_kkg_pi2_gridap_helper.jl` — voter (3), the Gridap FEM
  Laplace figure-helper; `figures/test_kkg_pi2_gridap.jl` — its
  figures-env verification.
- ADR-0001 (four-layer / additive architecture), ADR-0002 (BigFloat
  linear algebra via `GenericLinearAlgebra`).
- KKG 2015 §7 (Kapaev–Kitaev–… tritronquée figures 7.4/7.5) — the
  sector-fill-by-Laplace-solve precedent citing Trefethen.
- CLAUDE.md Rules 1, 2, 4, 5, 6, 9, 10.
