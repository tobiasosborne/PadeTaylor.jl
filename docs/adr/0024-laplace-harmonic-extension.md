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

  3. **Gridap.jl FEM Laplace solve** — a `PadeTaylorGridapExt` weak-dep
     package extension (ADR-0003 extensions pattern; separate bead
     `padetaylor-0ln.36`). An unstructured-FEM solution of the same
     Dirichlet Laplace problem, an algorithmically disjoint third voter.

The three interior fields are reconciled by a **per-grid-point majority
vote / agreement map** (assembled in F3, bead `padetaylor-0ln.33`):
points where the three methods agree to tolerance are accepted; the
agreement map is itself a figure-quality diagnostic. This mirrors the
project's standing triple-oracle discipline (CLAUDE.md Rule 4; cf. the
`SharedPade` Padé routine cross-validated against Chebfun's
`padeapprox.m`, a closed-form oracle, and a finite-difference oracle).
One numerical method can hide a systematic bug behind a plausible-looking
surface; three independent methods that agree cannot all share it.

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

- Bead `padetaylor-0ln.36` adds the `PadeTaylorGridapExt` weak-dep
  extension (third voter) with a `Project.toml` `[weakdeps]` line and a
  justification per ADR-0003 + Law 2.

- F3 (bead `padetaylor-0ln.33`) composes the three methods: it applies
  the `w = log x` conformal map, supplies the tritronquée boundary data,
  runs all three solvers, and assembles the per-grid-point majority vote
  and the agreement map.

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
- ADR-0003 — the weak-dep extension pattern that bead `0ln.36`'s Gridap
  voter (3) follows.
- ADR-0001 (four-layer / additive architecture), ADR-0002 (BigFloat
  linear algebra via `GenericLinearAlgebra`).
- KKG 2015 §7 (Kapaev–Kitaev–… tritronquée figures 7.4/7.5) — the
  sector-fill-by-Laplace-solve precedent citing Trefethen.
- CLAUDE.md Rules 1, 2, 4, 5, 6, 9, 10.
