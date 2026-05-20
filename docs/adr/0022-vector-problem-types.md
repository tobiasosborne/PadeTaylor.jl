# ADR-0022 — Vector problem types: the `NoumiYamadaProblem` and `PainleveHierarchyProblem` builders

**Status**: Accepted (2026-05-20) | **Beads**: `padetaylor-0ln.12` (V5a —
`NoumiYamadaProblem` builder + `A_n^(1)` RHS), `padetaylor-0ln.15` (V6 —
`PainleveHierarchyProblem` + `P_I^(2)` 4-vector RHS; sub-decision 4) |
**Plan**: `docs/v0p2_plan.md` rows V5a, V6.

## Decision

`NoumiYamada` — a new module — is the first of the v0.2 **vector problem
types**: builder layers that assemble a `VectorProblems.VectorPadeTaylorProblem`
for a specific *family* of vector ODEs, the way the scalar `Painleve`
module assembles a `PadeTaylorProblem` for the six Painlevé equations.
It exports two things:

  - `noumi_yamada_rhs(n; α)` — the right-hand-side factory for the
    even-parity Noumi–Yamada system `A_{2n}^{(1)}`; returns the closure
    `f(t, fvec) -> Vector` of length `2n+1`.
  - `NoumiYamadaProblem(n; α, f0, tspan, order)` — the builder: a small
    wrapper struct carrying the `(n, α)` system metadata plus the
    underlying `VectorPadeTaylorProblem`, with a forwarding
    `vector_solve_pade` method.

Three sub-decisions are recorded below.

### 1 — Per-system builder, not a value-function

`NoumiYamada` adds **no** `noumi_yamada(...)` value-function that would
return "the" solution.  A Noumi–Yamada system, like a Painlevé
equation, has infinitely many solutions; "which solution" must stay an
explicit caller choice expressed through the initial condition `f0`
(CLAUDE.md Rule 1, and ADR-0006's identical reasoning for `Painleve`).
What the module adds is a single discoverable *constructor* that selects
the correct `A_{2n}^{(1)}` RHS and assembles the vector IVP.

The builder is a **wrapper struct** (`NoumiYamadaProblem`, fields
`n`, `α`, `problem`) rather than a bare `VectorPadeTaylorProblem`.  This
mirrors `PainleveProblem` exactly: the wrapper makes the problem
self-describing (the system index `n` and parameters `α` are queryable
without re-deriving them from the closure), and a thin forwarding
`vector_solve_pade(prob::NoumiYamadaProblem; …)` method unwraps to
`prob.problem` — the same forwarding shape `Painleve` uses for
`solve_pade`.  Unlike the branch-point Painlevé equations, a
Noumi–Yamada system is integrated in its **natural `t`-frame**, so the
forwarding method carries no coordinate transform — it is a pure
unwrap-and-delegate.

### 2 — Full `(2n+1)`-component integration; constraint as a checked invariant

The even-parity system (NY1998 `\eqref{A2n}`,
`references/tex/noumi_yamada/NoumiYamada1998_higher_painleve_A1l_FunkEkv41/main.tex:85-88`)
is

    A_{2n}^{(1)}:   f_j' = f_j ( Σ_{1≤r≤n} f_{j+2r-1}
                                 −  Σ_{1≤r≤n} f_{j+2r} ) + α_j ,

a `(2n+1)`-component first-order vector ODE with parameters
`α_0, …, α_{2n}` summing to `1` (the `k = 1` normalisation; pillar B
`docs/v0p2_pillarB_noumi_yamada_findings.md` §1.2).  For that
normalisation the dependent variables satisfy `Σ_j f_j = t`.

We integrate the **full `(2n+1)`-component system directly** — no
constraint elimination, no reduction to the `2n` canonical Hamiltonian
pairs `(q_i, p_i)`.  The constraint `Σf_j = t` is an **invariant of the
flow**, not an externally imposed algebraic relation: summing the
system, every product `f_a f_b` in `Σ_j f_j·bracket_j` appears once with
`+` and once with `−` across the cyclic family and cancels, leaving
`Σ_j f_j' = Σ_j α_j = 1`, hence `d/dt (Σf_j − t) = 0`.

The constraint is therefore **checked, not eliminated**, in two places:

  - `NoumiYamadaProblem` *validates* `Σf0 = tspan[1]` at construction,
    so the conserved quantity `Σf_j − t` reads as `Σf_j = t` (not
    `Σf_j = t + const`) along the whole trajectory — it throws if the
    IC violates this (Rule 1);
  - it is the natural numerical self-test of a solved trajectory (test
    NY.1.4 asserts `Σ_j f_j(t) = t` at every breakpoint and at a dense
    interior point).

**Why full integration over canonical-pair elimination.** Keeping all
`2n+1` components makes the RHS a *verbatim* transcription of
`\eqref{A2n}` — the honest representation, with the cyclic structure
visibly intact and one obvious self-test (the constraint invariant).
Eliminating the constraint to integrate the `2n`-dimensional reduced
system would buy one fewer component at the cost of a bespoke
coordinate change per `n`, a non-obvious RHS, and the loss of the
constraint as a free correctness check.  The canonical-pair reduction
is the V6 hierarchy layer's concern, where the Hamiltonian structure is
the point; V5a stays close to the source equation.

### 3 — Even parity only for v0.2

`NoumiYamada` implements **only** the even-parity `A_{2n}^{(1)}` system.
This covers every v0.2 acceptance item: `A_2^{(1)}` (`n = 1`, the V5b
scalar-PIV backward-compat target) and `A_4^{(1)}` (`n = 2`, the V8a
pole-field figure and the PRD `n ≥ 4` acceptance).

The odd-parity `A_{2n+1}^{(1)}` system (NY1998 `\eqref{A2n+1}`,
`main.tex:92-102`) is the *quadratic-form* variant — its bracket is a
double sum of products `f_a f_b`, plus a linear `(k/2 − Σα)·f_j` term —
and reduces to PV (`A_3`) / a 4th-order system (`A_5`).  No v0.2
acceptance item needs it, so it is deliberately deferred.  Per Rule 9
the deferred corner is registered as bead `padetaylor-tji` with the
exact forcing condition: *needed if a v0.2+ target requires a PV-family
`A_{2n+1}` system*.  It would be implemented as an odd-parity branch of
this same module.

### 4 — `PainleveHierarchyProblem` — the second vector problem type (V6)

The V6 bead `padetaylor-0ln.15` adds `PainleveHierarchy`, the second
vector problem type, extending the pattern recorded above to the
**Painlevé-I hierarchy**.  It exports:

  - `painleve_hierarchy(:I, m; t)` — the companion-form RHS factory for
    the `m`-th PI-hierarchy member.  `m = 1` returns the 2-vector PI
    closure `y' = (y_2, 6 y_1² + x)`; `m = 2` returns the **4-vector
    `P_I^(2)` closure** `y' = (y_2, y_3, y_4, −10 y_2² − 20 y_1 y_3 −
    40(y_1³ − 6t y_1 + 6x))` — the genuine vector (4-component)
    Painlevé-type system and a v0.2 north-star target (KKG 2015
    eq. (p12), verified verbatim against the TeX source).
  - `PainleveHierarchyProblem(m; y0, xspan, t, order)` — the wrapper
    builder, carrying the `(m, t)` family metadata plus the underlying
    `VectorPadeTaylorProblem`, with a forwarding `vector_solve_pade`
    method.  Validates `length(y0) == 2m` and a non-degenerate `xspan`.
  - `pI2_tritronquee_ic(x0; t, n_terms)` — the leading-order tritronquée
    asymptotic-series seed for `P_I^(2)`.

**Companion-form vectorisation of a high-order scalar ODE.** This is the
structural distinction from `NoumiYamada`: a Noumi–Yamada system is
*already* a first-order vector ODE, whereas a PI-hierarchy member
`P_I^(m)` is a scalar ODE of order `2m`.  `painleve_hierarchy` reduces it
to first order by the standard companion form — state vector `y =
(u, u', …, u^{(2m-1)})`, the trivial chain `y_k' = y_{k+1}` for the first
`2m−1` rows, and the ODE solved for the top derivative in the last row.
The independent variable `x` is what the integrator steps in; the
parameter `t` (which enters `P_I^(2)` parametrically) is *captured into
the RHS closure*, not carried as a state component.  Otherwise the shape
is identical to sub-decision 1: a per-family builder producing a
`VectorPadeTaylorProblem`, a self-describing wrapper struct, and a
forwarding `vector_solve_pade` method (a new dispatch, not an edit).

**Scope — m ∈ {1, 2} explicit; m ≥ 3 deferred.** Members `m = 1` (PI, the
v0.1 reduction/consistency anchor — the `m = 1` companion RHS uses the
*same* `6u² + x` normalisation as v0.1 `PainleveProblem(:I)`, so no
rescaling is needed and the trajectories coincide) and `m = 2`
(`P_I^(2)`, the primary deliverable) are implemented **explicitly**, as
verbatim transcriptions of the source equations.  The general-`m` member
would require the Lenard–Gelfand–Dikii recursion to generate the
order-`2m` ODE symbolically — not a v0.2 acceptance item.
`painleve_hierarchy(:I, m)` throws for `m ≥ 3`; per Rule 9 the deferred
corner is registered as bead `padetaylor-qi0` with the forcing
condition recorded (*needed if a v0.2+ target requires a PI-hierarchy
member of order ≥ 6*).  This mirrors sub-decision 3's even-parity-only
scoping for `NoumiYamada` exactly.

## Context

v0.2 lifts PadeTaylor from scalar Painlevé equations to vector /
multi-component systems.  The lower layers are in place — `SharedPade`
(shared-`Q` Hermite–Padé, ADR-0019), `VectorCoefficients` (vector
Taylor jet, ADR-0020), `VectorStepper`, `VectorStepControl` (ADR-0021),
and `VectorProblems` (the `vector_solve_pade` driver).  What was missing
was a *problem-type* layer: a discoverable, self-describing constructor
for the specific vector systems the v0.2 north-star targets — the
Noumi–Yamada `A_n^{(1)}` higher-order Painlevé systems and (V6) the
`P_I^{(2)}` Painlevé-hierarchy companion form.

`NoumiYamada` is the **first** of these new vector problem types.  The
V6 Painlevé-hierarchy problem type (`painleve_hierarchy(:I, m)` /
`P_I^{(2)}`, bead `padetaylor-0ln.15`) will extend the same pattern: a
per-family builder producing a `VectorPadeTaylorProblem`, a wrapper
struct carrying the family metadata, and a forwarding
`vector_solve_pade` method.  Recording the builder shape here, before
V6 builds on it, keeps the two vector problem types consistent.

## Alternatives considered

  - **Return a bare `VectorPadeTaylorProblem`** (no wrapper struct).
    Rejected: it discards the `(n, α)` metadata, so a downstream
    consumer (the V5b PIV-reduction self-test, the V8a figure) would
    have to re-derive the system index from the closure.  The wrapper
    is one small struct and matches `PainleveProblem` — consistency
    across the scalar and vector builder layers is worth it.
  - **Eliminate the `Σf_j = t` constraint, integrate the reduced
    `2n`-dimensional system.** Rejected for V5a — see sub-decision 2.
  - **Implement both parities now.** Rejected — see sub-decision 3; no
    v0.2 acceptance item needs the odd case, and Rule 9 says document
    the deferred corner with its forcing condition rather than
    gold-plate.

## Consequences

- `NoumiYamada.jl` is a new module (87 effective LOC, well under the
  Rule-6 200-LOC cap; literate top docstring transcribing `\eqref{A2n}`,
  the cyclic structure, and the constraint-as-invariant).  It composes
  `VectorProblems` and loads after it.  Additive — no v0.1 *or* prior
  v0.2 module behaviour changes; `vector_solve_pade` gains one
  forwarding method (a new dispatch, not an edit).
- The vector RHS contract (ADR-0020) is honoured: `noumi_yamada_rhs`
  returns a closure `f(t, fvec)` of plain array + scalar arithmetic, so
  it composes with `VectorCoefficients` unchanged; `fvec` may carry
  `Taylor1` coefficients.
- Generic in `T`: the RHS and builder are element-type-agnostic; the
  `Σα = 1` and `Σf0 = t0` tolerances scale with `eps(float(real(T)))`,
  so a `BigFloat` parameter set is held to a tighter standard.
- Fail-loud (Rule 1): `n < 1`, an `α` / `f0` of wrong length, an `α`
  not summing to `1`, a state vector of wrong length passed to the
  RHS, and an `f0` violating `Σf0 = tspan[1]` all throw `ArgumentError`
  with a `Suggestion` line.
- The odd-parity `A_{2n+1}^{(1)}` system is a registered deferred bead
  (`padetaylor-tji`), not a silent gap.

## References

- `src/NoumiYamada.jl` — the V5a implementation (literate, 87 LOC).
- `test/noumi_yamada_test.jl` — the V5a suite (102 assertions, GREEN)
  + the NY.1.7 mutation-proof record (M1/M2/M3, all bit).
- `src/PainleveHierarchy.jl` — the V6 implementation (literate, 79 LOC):
  `painleve_hierarchy(:I, m)`, `PainleveHierarchyProblem`,
  `pI2_tritronquee_ic` (sub-decision 4).
- `test/painleve_hierarchy_test.jl` — the V6 suite (59 assertions,
  GREEN) + the PH.1.6 mutation-proof record (M1/M2/M3, all bit).
- `references/tex/painleve_hierarchy/KapaevKleinGrava2015_PI2_tritronquee_ConstrApprox41/tritronquee_coeff.tex`
  lines 124-130 — PI and the `P_I^(2)` eq. (p12), the verbatim source of
  the V6 companion RHSs.
- `docs/v0p2_pillarC_painleve_hierarchy_findings.md` — §1 (the explicit
  `P_I^(2)` equation), §2 (the 4-vector companion form), §3 (the PI
  hierarchy / Lenard recursion), §4 (the tritronquée asymptotic IC).
- Bead `padetaylor-qi0` — the deferred general-`m` (m ≥ 3) PI-hierarchy
  member via the Lenard–Gelfand–Dikii recursion.
- `references/tex/noumi_yamada/NoumiYamada1998_higher_painleve_A1l_FunkEkv41/main.tex`
  lines 85–88 — equation `\eqref{A2n}`, the verbatim source of the
  even-parity RHS; lines 92–102 — `\eqref{A2n+1}`, the deferred
  odd-parity system.
- `docs/v0p2_pillarB_noumi_yamada_findings.md` — §1 (the system, both
  parities, the dimension table), §1.2 (the `Σf_j = t` constraint),
  §1.4 (small-`l` explicit formulas), §5.4 / §7.2 (the `A_4` Type C
  rational oracle `f_j = t/5`).
- `src/VectorProblems.jl` — `VectorPadeTaylorProblem` /
  `vector_solve_pade`, the driver this builder feeds.
- `src/Painleve.jl` — the scalar per-equation builder pattern mirrored
  here (ADR-0006).
- ADR-0019 — the shared-`Q` Padé; ADR-0020 — the vector Taylor jet and
  the RHS contract; ADR-0021 — the vector step control.
- `docs/v0p2_plan.md` — §"The 10 phases", V5a row; V6 row (the next
  vector problem type, which extends this pattern).
- CLAUDE.md Rules 1, 4, 6, 9, 10.
