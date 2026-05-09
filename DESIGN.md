# DESIGN.md — PadeTaylor.jl Stage 1 design (Julia v1)

> **Status: Stage 1 deliverable.** Mandates the Julia-side architecture
> for v1; the TS-mirror design (`@workbench/pade-taylor`) is deliberately
> deferred to Stage 3 ("get Julia working 100%, *then* spec the
> workbench version"). Every claim cites either `RESEARCH.md §X.Y` or
> a primary-literature line in `references/markdown/<file>.md`.

## 0. Discipline (read every session)

Inheriting verbatim from scientist-workbench's two laws:

1. **Law 1 — Ground truth before code.** No code without an open ADR /
   open paper / open existing-file tab. Citations are by file path
   (e.g. `references/markdown/FW2011_painleve_methodology_JCP230/
   FW2011_painleve_methodology_JCP230.md:128`), not paraphrase from
   memory.
2. **Law 2 — Docs in lockstep.** Every algorithmic decision has an
   ADR under `docs/adr/NNNN-title.md`; every module has a
   docstring-as-chapter at the top of its source file.

**TDD discipline (workbench Rule 6).** Two valid shapes:

- **Port-and-verify** — apply to `RobustPade` (port of GGT 2013
  Algorithm 2 + Chebfun reweighting), `StepControl` (port of Jorba-Zou
  2005 §3.2), and the `coefficient-via-recurrence` paths. Procedure:
  port faithfully → write invariant tests → **mutation-prove** (perturb
  impl, confirm RED, restore) → cross-validate against an oracle when
  one exists.
- **Spec-from-scratch** — apply to `PadeStepper` (the four-layer
  pipeline integrator), the `Problems` IVP shape, and the
  `step_complex` primitive. Procedure: classic RED → GREEN → commit.

**Per-file ≤ ~200 LOC** (workbench Rule 12 inheritance). When a module
exceeds this, split.

**Beads** is the only tracker. `bd create` for any friction, problem,
deferral, or follow-up. `bd close <ids>` at session end.

**No parallel Julia agents** — precompile-cache conflicts. Single-
session Julia work; subagents OK only for read-only literature work.

**Fail fast, fail loud.** No silent fallbacks. `ToolError`-style
exceptions with `suggestion` text, not "return zero on error."

## 1. Module decomposition

Six modules, dependency order leaf-up:

```
                LinAlg            ← Tier 0 (leaf)
                   │
        ┌──────────┼─────────┐
        │          │          │
    RobustPade  Coefficients (uses TaylorSeries.jl)
                                  ← Tier 1
                   │
              StepControl   ← Tier 2 (uses Coefficients output type
                              + Polynomials.jl::roots for FW heuristic)
                   │
              PadeStepper   ← Tier 3 (combines all four)
                   │
               Problems     ← Tier 4 (public API)
                   │
        CommonSolveAdapter  ← Tier 5 (optional SciML integration)
```

| Module | LOC budget | Purpose |
|---|---|---|
| `LinAlg` | 60 | Single dispatch: `pade_svd(A::Matrix{T})` → `(U, S, Vt)` for `T ∈ {Float64, BigFloat, Arblib.Arb}`. Routes to `LinearAlgebra.svd` for `Float64`; `GenericLinearAlgebra.svd` over `Matrix{BigFloat}` otherwise. *Why a wrapper at all?* (a) Single point to add Arb-conversion shim; (b) stable internal API that the Padé module can rely on without caring about Julia stdlib evolution; (c) the place to add the optional fallback to a hand-rolled Demmel-Veselić Jacobi if `GenericLinearAlgebra` proves flaky on near-deficient GGT matrices. |
| `RobustPade` | 180 | GGT 2013 Algorithm 2 + Chebfun's QR-reweighting. `robust_pade(c, m, n; tol)` → `PadeApproximant{T}` containing `a, b, μ, ν` and rescaling info. **The single point of truth for our Padé conversion.** |
| `Coefficients` | 80 | Wrap `TaylorSeries.jl::Taylor1{T}` for our use. `taylor_coefficients(f, z0, y0, order)` → `Vector{T}` for first-order ODE. `taylor_coefficients_2nd(f, z0, y0, y1, order)` → `Vector{T}` for second-order (PI/PII/PIV are all 2nd-order). |
| `StepControl` | 150 | Two strategies: `step_jorba_zou(coefs, ε, ε_r)` (default; per Jorba-Zou 2005 §3.2 eq. 3-8) and `step_pade_root(P, Q, z_current, target)` (FW 2011 §3.1 path heuristic, refined by denominator-root distance). Strategy enum: `:jorba_zou`, `:pade_root`. |
| `PadeStepper` | 200 | The four-layer pipeline: Taylor coefs → Padé conversion → step direction & size → take step. `pade_step(state, problem)` → `(z_new, y_new)`. Handles direction selection (5-direction sampling vs steepest-descent). |
| `Problems` | 150 | `PadeTaylorProblem{F, T}(f, y0, zspan; order=30)` and `solve_pade(prob; kwargs...)` → `PadeTaylorSolution`. Trajectory storage; dense interpolation via stored Padé approximants. |
| `CommonSolveAdapter` | 50 | Optional. `CommonSolve.solve(prob, alg; kwargs...) → solve_pade(...)`. Import only when `CommonSolve.jl` is in user's Project.toml; non-essential dep. |

Total v1 source: **~870 LOC**. Tests + docstrings: **~700 LOC**.
Total v1: **~1570 LOC**. Tractable in 2–4 weeks of focused work.

**Public API surface** (exported from `PadeTaylor` module):
```
PadeTaylorProblem, solve_pade, PadeTaylorSolution
robust_pade, PadeApproximant         # for advanced users
taylor_coefficients
```

Internal modules (`PadeTaylor.LinAlg`, etc.) are intentionally not
re-exported; users go through the public API.

## 2. Dependencies

```toml
[deps]
TaylorSeries          = "..."  # validated empirically per RESEARCH.md §3.3
GenericLinearAlgebra  = "..."  # bigfloat SVD (per §5.1, §5.3)
Polynomials           = "..."  # roots of Q for step_pade_root

[weakdeps]
Arblib                = "..."  # Arb element type (extension package)
CommonSolve           = "..."  # SciML adapter (extension package)
```

**Extensions pattern (Julia 1.9+):**

- `PadeTaylorArblibExt.jl` — pulled in only when user imports `Arblib`.
  Shims `Arblib.Arb` → `BigFloat` conversion for `pade_svd`, plus the
  `radius` reporting in `PadeTaylorSolution` for honest-error-bar mode.
- `PadeTaylorCommonSolveExt.jl` — pulled in only when user imports
  `CommonSolve`. Implements the `solve(prob, alg; kwargs...)` adapter.

Justification: keeps the core dep graph minimal; users who want Arb
or CommonSolve opt in. Pattern is canonical in modern Julia packages
(SciML, JuliaDiff, etc.).

## 3. Test corpus — what every module must answer

Three target tiers, each with an explicit cross-validation oracle.

### 3.1 Tier A: Float64 baseline (load-bearing for v1)

Source: FW 2011 Table 5.1 verbatim
(`references/markdown/FW2011_painleve_methodology_JCP230/
FW2011_painleve_methodology_JCP230.md:385–391`).

Test problem: Weierstrass ℘ — `u'' = 6u²` with `c₁ = -1, c₂ = 2`,
giving `u(0) ≈ 1.071822516416917, u'(0) ≈ 1.710337353176786`.

Reference values:
- `u(30) = 1.095098255959744`
- `u(10⁴) = 21.02530339471055`
- `u(28.261) = 9.876953517025014 × 10⁶`

Acceptance: `(order, h) = (30, 0.5)`, relative errors:
- at `z = 30`: `|err| ≤ 1e-13` (FW achieve `7.62e-14`; we accept up to
  ~2× their result).
- at `z = 10⁴`: `|err| ≤ 5e-10` (FW: `2.34e-10`).
- at `z = 28.261`: `|err| ≤ 2e-9` (FW: `7.92e-10`).

### 3.2 Tier B: BigFloat / arb-prec (load-bearing for arb-prec correctness)

Same Weierstrass ℘ test problem. Reference values computed at much
higher precision via the closed-form ℘-function (Wolfram
`WeierstrassP[z, {0, 2}]` gives 100+ digits).

At precision 256 bits (~77 decimal digits):
- `u(30)`: relative error ≤ `1e-65` (well above the `2⁻²⁰⁰` baseline
  the empirical probe demonstrated for raw coefficient arithmetic).

### 3.3 Tier C: PI tritronquée IC (verification of "general purpose")

Source: FW 2011 eq. (4.1) `ref/markdown/FW2011_*.md:226`.

`u(0) ≈ -0.1875543083404949, u'(0) ≈ 0.3049055602612289` integrating
PI `u'' = 6u² + z`. Acceptance: produce a pole-field figure that
qualitatively reproduces FW Fig. 3.1 — sectorial structure, pole
density `O(R^{5/2})` within radius `R`, near-tritronquée appearance.

### 3.4 Tier D: Pure GGT 2013 Padé routine validation

Source: GGT 2013 Figs. 5–9 + the `padeapprox.m` reference impl.

Test cases (parallel run: our `robust_pade` vs MATLAB-`padeapprox.m`
output captured ahead of time):
- `f(z) = exp(z)`, type `(20, 20)`. Standard easy case.
- `f(z) = tan(z⁴)`, type `(100, 100)`. Froissart-doublet stress test.
- `f(z) = log(1.2 - z)`, type `(20, 20)`. Branch-cut behaviour.
- `f(z) = 1 + z²`, type `(1, 1)`. Defect-1 ill-posedness sanity.

For each: our `robust_pade(c, m, n)` must produce `(μ, ν)` matching
Chebfun's, with `‖a - a_chebfun‖_∞ ≤ tol_match` and `‖b - b_chebfun‖_∞
≤ tol_match`. Capture Chebfun's outputs via Octave/MATLAB if available
locally; if not, transcribe from GGT 2013 figures and accept slightly
weaker bounds.

## 4. Granular execution plan — module-by-module, leaf-up

Each step is a discrete commit. **Do not start the next step until
the current step is GREEN with mutation-proven tests.**

### Phase Z — pre-conditions (no code yet)

1. **ADR-0001: four-layer architecture.** Document the
   Coefficients-Padé-Step-Solver decomposition. Cite FW 2011 §1–§5
   and GGT 2013 §2 for each layer.
2. **ADR-0002: bigfloat-SVD via GenericLinearAlgebra.** Per
   `RESEARCH.md §5.1` — argument for one-sided Jacobi over
   Demmel-Kahan, plus the `Arb → BigFloat` conversion shim and its
   precision-loss caveat.
3. **ADR-0003: extensions pattern.** When to use Pkg.jl extensions
   vs. hard deps. Argument: keep core deps minimal (TaylorSeries +
   GenericLinearAlgebra + Polynomials), gate Arb and CommonSolve
   behind extensions.
4. **Replace generic CLAUDE.md** with PadeTaylor-specific rules
   (workbench Rule 9 inheritance: ground-truth, ≤200 LOC,
   no parallel Julia agents, beads-only, fail-fast).
5. **`Project.toml` skeleton + `src/PadeTaylor.jl` umbrella.**
   Empty submodules, no implementations yet. Confirm `julia
   --project=. -e 'using Pkg; Pkg.instantiate(); using PadeTaylor'`
   loads cleanly.

### Phase 1 — `LinAlg` (Tier 0, ~60 LOC + ~80 LOC tests)

**Dependencies**: `LinearAlgebra` (stdlib), `GenericLinearAlgebra`.

**Public API**:
```julia
"""
    pade_svd(A::AbstractMatrix{T}) -> (U, S, Vt)

SVD with relative accuracy on small singular values, regardless of κ(A).

For `T = Float64`: dispatches to LAPACK Demmel-Kahan via
`LinearAlgebra.svd` (GGT 2013's tolerance regime tolerates this at
Float64; see `RESEARCH.md §5.3`).

For `T <: AbstractFloat` other than `Float64` (incl. `BigFloat`):
dispatches to `GenericLinearAlgebra.svd` (one-sided Jacobi
Demmel-Veselić; relative-accuracy guarantee `c · 2⁻ᵖ · σᵢ` per SV;
load-bearing for arb-prec rank-counting in `RobustPade`).
"""
function pade_svd end
```

**TDD step 1.1**: write `test/linalg_test.jl` with these RED-then-GREEN
tests:

| test | input | expected | shape |
|---|---|---|---|
| 1.1.1 | `Matrix{Float64}` random `5 × 6`, seed 42 | matches `LinearAlgebra.svd(A)` exactly | spec-from-scratch |
| 1.1.2 | `Matrix{Float64}` Hilbert-10 | reconstructs `A` to `1e-9` rel error (well-conditioned regime) | spec-from-scratch |
| 1.1.3 | `Matrix{BigFloat}` Hilbert-10 at `precision=256` | reconstructs to `1e-65` rel error | port-and-verify (against `GenericLinearAlgebra.svd`) |
| 1.1.4 | `Matrix{BigFloat}` rank-deficient `[1 2; 2 4]` | `S[2] / S[1] < 2⁻¹⁰⁰` (genuine zero) | port-and-verify |
| 1.1.5 | mutation-proof: invert a 2×2 SVD's `Vt` row → assert test 1.1.3 fails | RED | mutation-proof |

**Acceptance for Phase 1**: all 5 tests pass; `bd close
<linalg-bead>`; commit message references which RESEARCH.md sections
the design draws from.

### Phase 2 — `RobustPade` (Tier 1, ~180 LOC + ~250 LOC tests)

**Dependencies**: `LinAlg`.

**Source-of-truth for the port**: `external/chebfun/padeapprox.m`
(commit `7574c77`, lines 263–296). Mirror **the algorithm**, not the
MATLAB idioms.

**Public API**:
```julia
"""
    PadeApproximant{T}

Output of `robust_pade`. Fields: `a::Vector{T}` (numerator coefs),
`b::Vector{T}` (denominator coefs, with `b[1] = 1`), `μ::Int`
(exact numerator degree), `ν::Int` (exact denominator degree),
`γ::T` (rescaling parameter; default 1).
"""
struct PadeApproximant{T} ... end

"""
    robust_pade(c::AbstractVector{T}, m::Int, n::Int; tol = default_tol(T))
        -> PadeApproximant{T}

Robust Padé approximant of type (m, n) to the Taylor series with
coefficients `c[1] = c₀, c[2] = c₁, …, c[m+n+1] = c_{m+n}`.

Implements GGT 2013 Algorithm 2 (`references/markdown/GGT2013_*.md:
219-235`) with the Chebfun reweighting trick (lines 236-241).

`tol` default:
- `Float64`: `1e-14` (GGT default for Float64 rounding-error problems).
- `BigFloat`/`Arb`: `2.0^(-precision(T) + 10)` (10 bits of slack
  above the working precision floor; `RESEARCH.md §2.1.3`).

Throws `PadeError` if `c` is too short (`length(c) < m + n + 1`),
all `c[1..m+1]` are below `tol`, or any input is non-finite. Returns
the unique normalized minimal-degree representation per
`Theorem 2.1` of GGT 2013.
"""
function robust_pade end

"""
    evaluate(P::PadeApproximant{T}, z::Number) -> T

Evaluate `P` at `z` via Horner's rule on numerator and denominator.
Throws `DomainError` if denominator vanishes within `tol` of zero
at `z` (callers should catch and shorten step).
"""
function evaluate end
```

**TDD step 2.1**: write `test/robustpade_test.jl` with these:

| test | input | expected | shape |
|---|---|---|---|
| 2.1.1 | `c = [1, 1, 1/2, 1/6, 1/24, 1/120]`, `(m, n) = (2, 2)` (exp) | `(μ, ν) = (2, 2)`; `b ≈ [1, -1/2, 1/12]`; `a ≈ [1, 1/2, 1/12]` | spec-from-scratch (closed-form Padé of exp) |
| 2.1.2 | `c = exp coeffs to order 40`, `(m, n) = (20, 20)` | matches Chebfun `padeapprox` output coefs to `1e-12` rel | port-and-verify (oracle: pre-computed Chebfun output) |
| 2.1.3 | `c = log(1.2 - z) coeffs to order 40`, `(m, n) = (20, 20)` | `(μ, ν)` matches Chebfun (degeneracy reduction triggers) | port-and-verify |
| 2.1.4 | `c = tan(z⁴) coeffs`, `(m, n) = (20, 20)` | poles of `b` cluster on the 8 expected rays; no Froissart doublets | port-and-verify (visual + count) |
| 2.1.5 | `c = [1, 0, 1]`, `(m, n) = (1, 1)` | defect-1 case; output `(μ, ν) = (0, 0)`, `r ≡ 1` (per GGT 2013 eq 7.1) | spec-from-scratch |
| 2.1.6 | `c = noisy [1, 1, 1, …]` with `1e-6` Gaussian noise; `tol = 1e-5` | recovers `(μ, ν) = (0, 1)` (`1/(1-z)`); GGT 2013 Fig. 4 right panel | port-and-verify |
| 2.1.7 | mutation-proof: comment out the QR-reweighting block → assert 2.1.4 fails | RED | mutation-proof |
| 2.1.8 | `T = BigFloat(precision=256)`; `c = exp coeffs to order 40`; `(m, n) = (20, 20)` | matches Float64 `b` to `1e-12`; coefs themselves accurate to `2⁻²⁰⁰` | port-and-verify (precision-tier) |

**Acceptance for Phase 2**: all 8 tests pass; `bd close
<robustpade-bead>`; commit message references GGT 2013 + `padeapprox.m`
line numbers.

### Phase 3 — `Coefficients` (Tier 1, ~80 LOC + ~120 LOC tests)

**Dependencies**: `TaylorSeries.jl`.

**Public API**:
```julia
"""
    taylor_coefficients_1st(f, z0::T, y0::T, order::Int) where T -> Vector{T}

For first-order ODE `dy/dz = f(z, y)`, compute the Taylor coefficients
`c₀, c₁, …, c_order` of `y(z₀ + h)` such that
`y(z₀ + h) ≈ c₀ + c₁ h + c₂ h² + … + c_order h^order`.

Implements method (b) from FW 2011 §2.1.2 (`references/markdown/
FW2011_*.md:96-107`): substitute the truncated Taylor expansion into
`y'(z₀+h) = f(z₀+h, y(z₀+h))` and bootstrap one coefficient per pass
via integration of the previous Taylor truncation. Mechanically: use
`TaylorSeries.jl::Taylor1{T}` operator overloading to compute
`f(z₀+h, y(h))` as a `Taylor1{T}` and integrate.

For second-order ODEs use `taylor_coefficients_2nd`.
"""
function taylor_coefficients_1st end

function taylor_coefficients_2nd end  # for u'' = f(z, u, u')
```

**TDD step 3.1**:

| test | input | expected | shape |
|---|---|---|---|
| 3.1.1 | `f(z, y) = y`, `z₀ = 0, y₀ = 1`, order 10 | `c[k] = 1/k!` (i.e. `exp(z)` Taylor) | spec-from-scratch |
| 3.1.2 | `f(z, y) = z² + y²`, `z₀ = 0, y₀ = 0`, order 14 (FW 2011 §2.2.1 demo) | matches `Taylor1` of `t * J_{3/4}(t²/2) / J_{-1/4}(t²/2)` to ball-equality at order 14 | port-and-verify (ref: FW 2011 §2.2.1) |
| 3.1.3 | `f(z, u, u') = 6u² + z` (PI), `z₀ = 0, u₀ = -0.1875, u'₀ = 0.3049`, order 30 | numeric coefs match `WeierstrassP`-derived expansion to `1e-12` | port-and-verify |
| 3.1.4 | `T = BigFloat(precision=256)`, same as 3.1.1, order 60 | radii `< 2⁻²⁰⁰` (per `RESEARCH.md §3.3` empirical) | port-and-verify (precision-tier) |
| 3.1.5 | mutation-proof: replace one rule in `TaylorSeries`'s `exp` with sign flip → assert 3.1.1 fails | RED | mutation-proof |

**Acceptance for Phase 3**: all 5 tests pass; `bd close <coefs-bead>`.

### Phase 4 — `StepControl` (Tier 2, ~150 LOC + ~120 LOC tests)

**Dependencies**: `Polynomials.jl::roots` for FW heuristic; nothing
else.

**Public API**:
```julia
"""
    step_jorba_zou(coefs::Vector{T}, eps_abs::Real; eps_rel = eps_abs)
        -> T

Estimated step length per Jorba-Zou 2005 §3.3.1 eq. 11, in the
fixed-order form ported verbatim from `TaylorIntegration.jl::stepsize`.
With p = `length(coefs) - 1`:
    h = min over k ∈ {p-1, p} of (ε / |c[k+1]|)^(1/k)
where `ε = eps_abs` if `eps_abs ≥ eps_rel · |c₀|` else `eps_rel · |c₀|`
(TI.jl ε-resolution).  Zero candidates are skipped; if both are zero,
fall back to the TI.jl `_second_stepsize` scan.

CORRECTED 2026-05-09 (commit fixing Phase 4): an earlier version of
this spec asserted `h = (ρ/e²)·exp(-0.7/(p-1))`.  The `0.7` constant
has no source in either the paper or TI.jl (grep both: zero matches);
it was a hallucination at spec-write time.  The paper's `/e²` is for
the asymptotically-optimal-(p,h) regime where p is also free; we fix
p = 30 (FW 2011), so the `/e²` is absorbed into the ε substitution
(see the worked algebra in `src/StepControl.jl`'s top-of-file
docstring).  Three-source consensus on the canonical test case
(c[k]=1/k!, p=30, ε=1e-12): TI.jl ≡ mpmath ≡ wolframscript =
4.50120637033898607690318848315848021108523516609 (47 decimal digits).
"""
function step_jorba_zou end

"""
    step_pade_root(P::PadeApproximant{T}, z_current::Number, target::Number)
        -> T

Step length to the nearest pole of `P` (root of denominator) along the
direction `(target - z_current)`. Implements FW 2011 §3.1's pole-
distance heuristic. Returns a step length less than the full
`|target - z_current|` if a pole intervenes.
"""
function step_pade_root end
```

**TDD step 4.1**:

| test | input | expected | shape |
|---|---|---|---|
| 4.1.1 | `coefs = [1/k! for k=0..30]`, `eps_abs = 1e-12` | `h = 4.501206370338986` (three-source consensus pinned in `test/_oracle_stepcontrol.jl`) | port-and-verify (paper §3.3.1 + TI.jl) |
| 4.1.2 | same input | round-trip equality with `TaylorIntegration.stepsize` to 1e-15 rel | port-and-verify |
| 4.1.3 | `P = 1/(1 - z/2)`, `z_current = 0`, `target = 5` | step = `2` (pole at `z = 2`) | spec-from-scratch |
| 4.1.4 | `P = 1/(1 + (z - 3)²)`, `z_current = 0`, `target = 5` | step = `3` (real-axis projection of `3 ± i` poles) | spec-from-scratch |
| 4.1.5 | mutation-proof: change exponent `1/k → 1/(k+1)` in step_jorba_zou → 4.1.1/4.1.2 RED; change `real → imag` in step_pade_root projection → 4.1.3/4.1.4 RED | RED | mutation-proof |

**Acceptance for Phase 4**: all 5 tests pass.

### Phase 5 — `PadeStepper` (Tier 3, ~200 LOC + ~200 LOC tests)

**Dependencies**: all of `LinAlg`, `RobustPade`, `Coefficients`,
`StepControl`. **This is where the whole pipeline first runs end-to-
end.**

**Public API (as shipped 2026-05-09)**:
```julia
mutable struct PadeStepperState{T}
    z::T
    u::T
    up::T
end

"""
    pade_step!(state::PadeStepperState{T}, f, order::Int, h::Real) -> state

Take one Padé-Taylor step of `u'' = f(z, u, u')`, mutating state in
place.  Step direction (5-direction / steepest-descent) and adaptive
step-size selection are deferred to Phase 6's solve_pade — Phase 5
accepts an explicit `h`.
"""
function pade_step! end
```

**Algorithm (5 steps, as shipped — corrected from DESIGN's original
6-step sketch)**:
1. Generate Taylor coefs of u via `Coefficients.taylor_coefficients_2nd`.
2. **Rescale `c̃_k = h^k · c_k`** (variable substitution `h' = h·t` per
   FW 2011 §3.2; absolutely necessary near poles where raw Taylor
   coefs span 30+ orders of magnitude and `RobustPade`'s
   `tol·‖c‖∞`-based zero-detection fires spuriously). Build diagonal
   Padé `P_u = robust_pade(c̃, m, n)` with `m = n = order ÷ 2`.
3. Evaluate `new_u = P_u(1)`.
4. **Compute `new_up` via analytic differentiation of P_u at t=1**,
   then divide by h (chain rule from `u(z+h') = P_u(h'/h)`).  NOT a
   second `robust_pade` call on the formal-derivative coefficients —
   empirically the re-Padé path loses ~40× accuracy on test 5.1.2.
5. Mutate state.

See `src/PadeStepper.jl`'s top docstring + `docs/worklog/003-phase-5-
padestepper.md` for the rationale on (a) the h^k rescaling and (b) the
analytic-differentiation choice.

**TDD step 5.1**:

| test | input | expected | shape |
|---|---|---|---|
| 5.1.1 | one ℘ step at `(z=0, h=0.5, order=30)` from FW ICs | `u(0.5), u'(0.5)` match ℘ closed-form to 1e-13 rel | port-and-verify (3-source oracle: closed-form ≡ NDSolve ≡ mpmath.odefun) |
| 5.1.2 | compose: continue 5.1.1 to z=0.9 (h=0.4) | matches ℘ closed-form to 1e-12 rel [DEVIATION: DESIGN said z=1.0, but z=1.0 IS the pole] | port-and-verify |
| 5.1.3 | one PI step at FW ICs (`u'' = 6u² + z`) | `u(0.5)` finite, matches NDSolve oracle to 1e-13 rel; differs from ℘ by ≈0.0296 | port-and-verify [DEVIATION: DESIGN said tritronquée IC; no paper-pinned tritronquée IC at high precision available, FW IC tests same pipeline + +z RHS] |
| 5.1.4 | near-pole step (z=0.9 → 0.95, h=0.05) on u''=6u² | matches closed-form to 1e-11 rel [DEVIATION: DESIGN said z=1.36 pole; the actual pole on this lattice is z=1, oracle confirms 1/u(1)=0] | port-and-verify |
| 5.1.5 | mutation-proof: (A) sign flip on Padé eval point bites all `state.u`; (C) chain-rule `*h_T` instead of `/h_T` bites all `state.up` | RED | mutation-proof |

**Acceptance for Phase 5**: all tests pass (16 individual @test calls
across 4 testsets, 207/207 total). **The inner-loop algorithm is
proven correct**; remaining phases are scaffolding.

### Phase 6 — `Problems` (Tier 4, ~150 LOC + ~150 LOC tests)

**Dependencies**: `PadeStepper`.

**Public API**:
```julia
"""
    PadeTaylorProblem{F, T}(f, y0, zspan; order=30, ...)

Definition of an analytic ODE IVP. `f` is a Julia function with
signature `f(z, y) -> dy/dz` (first-order) or `f(z, y, y') -> y''`
(second-order, autodetected from `y0`'s tuple-arity).

`y0`: initial condition. For a first-order ODE, `y0::T`; for second-
order, `y0::Tuple{T, T}` is `(u(z0), u'(z0))`.

`zspan::Tuple{T, T}`: integration domain endpoints. For complex
integration, both can be `Complex{T}`.
"""
struct PadeTaylorProblem ... end

"""
    solve_pade(prob::PadeTaylorProblem; step_strategy=:jorba_zou,
               h_max=0.5, tol=default_tol(T), max_steps=10_000)
        -> PadeTaylorSolution

Integrate `prob` from `zspan[1]` to `zspan[2]`. Stores trajectory
+ Padé approximants for dense interpolation. Throws `IntegrationError`
on `max_steps` exhaustion or persistent step-rejection.
"""
function solve_pade end

"""
    PadeTaylorSolution{T}

Result of `solve_pade`. Fields: `z::Vector{T}`, `y::Vector{Y}`,
`pade::Vector{PadeApproximant{T}}`. Callable: `sol(z) -> Y` for
dense interpolation via the stored Padé.
"""
struct PadeTaylorSolution ... end
```

**TDD step 6.1**:

| test | input | expected | shape |
|---|---|---|---|
| 6.1.1 | ℘ from `z=0` to `z=30`, FW 2011 ICs, `(order=30, h=0.5)` | `sol(30)` matches `1.095098255959744` to `5e-13` (FW Tier A) | port-and-verify (FW Table 5.1) |
| 6.1.2 | ℘ from `z=0` to `z=10⁴`, same | `sol(10⁴)` matches `21.02530339471055` to `5e-10` (FW: 2.34e-10) | port-and-verify |
| 6.1.3 | ℘ at `z = 28.261` (high on pole wall) | `sol(28.261) ≈ 9.876953517025014e6` to `2e-9` | port-and-verify |
| 6.1.4 | ℘ with `T = BigFloat(precision=256)` from `z=0` to `z=30` | `sol(30)` matches closed-form to `1e-65` | port-and-verify (Tier B) |
| 6.1.5 | dense interp: `sol(7.123)` for the 6.1.1 problem | matches closed-form ℘ to `5e-13` (Padé-interpolated, no extra step) | spec-from-scratch |
| 6.1.6 | mutation-proof: replace `pade_step!` body with no-op → 6.1.1 fails with `assert` mismatch | RED | mutation-proof |

**Acceptance for Phase 6**: all 6 tests pass; **at this point v1 is
algorithmically complete**.

### Phase 7 — `CommonSolveAdapter` extension (~50 LOC + ~30 LOC tests)

Optional. Implements the SciML `CommonSolve.solve` interface as a
thin wrapper. Acceptance: a SciML user calling `solve(prob,
PadeTaylorAlg(); kwargs...)` gets bit-identical output to direct
`solve_pade(prob; kwargs...)` for the Tier A test cases.

### Phase 8 — `PadeTaylorArblibExt` extension (~80 LOC + ~80 LOC tests)

Implements:
- `pade_svd(::Matrix{Arblib.Arb})` — converts via `Matrix{BigFloat}`
  (lossy: rounds Arb mid to BigFloat, discards radius). The SVD is
  then in the `BigFloat` regime.
- `default_tol(::Type{Arb})` — `Arblib.precision_bits` aware.
- `radius` reporting in `PadeTaylorSolution{Arb}` — passes through
  Arb's per-coefficient radius from `taylor_coefficients`.

**Acceptance**: all Phase-1–6 tests run with `T = Arb` substituted
where applicable.

### Phase 9 — Tier C qualitative reproduction (no new code)

Run PI tritronquée at FW 2011 eq. (4.1) ICs, generate pole-field
plot, compare visually to FW Fig. 3.1. Acceptance: sectorial structure
matches; pole density matches `O(R^{5/2})` count. **Closes v1.**

## 5. What we won't do in v1 (explicit non-goals — deferred to v2)

- Real-axis-network integration with the FW 2011 §3.1 branching path
  algorithm (`step_complex` is the v1 primitive; v2 builds the network
  on top).
- BVP coupling for smooth regions (FW 2011 §3.2). v2.
- Pole-field-edge detection (FW 2011 §3.2.2). v2.
- Hermite-Padé. v2.
- Riemann-surface tracking (FFW 2017). v2.
- Adaptive order. v2.
- Verified-arithmetic Cauchy-Hadamard truncation bounds. v1 flag /
  early v2.
- The TS-mirror `@workbench/pade-taylor`. **Stage 3** — only after
  Julia v1 is 100% working and the cross-validation oracle exists.

## 6. Cross-validation discipline (carries into Stage 3)

When the TS mirror lands (Stage 3), each test in §3.1–§3.4 of this
document gets a TS counterpart that produces bit-byte-identical output
on Float64 (per workbench's symbolic-tier determinism contract,
modulo platform fingerprint for numerical operations) and matches to
`N` digits at the chosen arb-prec setting.

**Cross-validation is the load-bearing reason both impls exist** —
neither alone is "the reference"; the *pair* is.

## 7. Open implementation-time questions (file beads as they arise)

From `RESEARCH.md §8`, repeated here for visibility:

1. GGT 2013 vs FW 2011 Toeplitz-direct empirical comparison — runs
   alongside Phase 6.
2. Auto-rescaling parameter `γ` per Fornberg 1981 — likely v1.5.
3. `Polynomials.jl::roots` accuracy with `Arb` coefs — relevant in
   Phase 4, may force a `GenericSchur.jl` companion-matrix path.
4. `GenericLinearAlgebra.svd` direct `Arb` element type — Phase 8;
   fall back to `BigFloat` conversion if it doesn't compile.
5. Symbolic AD path for coefficient generation — v2.
6. Step-direction selection — Phase 5 sub-decision.
7. Path-network self-consistency at API level — v2.

Each gets a beads issue when encountered.
