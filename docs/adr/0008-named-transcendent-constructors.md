# ADR-0008 — Named Painlevé transcendent constructors

**Status**: Accepted (2026-05-14) | **Bead**: `padetaylor-c86`
**Context**: Continues the Painlevé-infrastructure roadmap (ADR-0006
builder, ADR-0007 solution wrapper).  ADR-0006 explicitly deferred
"named special-solution constructors (`tritronquee(:I)`,
`hastings_mcleod(:II)`, …) — a convenience on top of this layer once it
exists".  ADR-0007 shipped that layer (`PainleveSolution`) with a
`name::Union{Symbol,Nothing}` field reserved for exactly this.  This
ADR cashes the deferral in.

## Context

A user reaching for "the Painlevé-I tritronquée" or "the
Hastings–McLeod solution" should not have to know the 16-digit initial
conditions that pin them.  Those numbers are the *definition* of the
named solution — they belong in the package, cited to the paper that
computed them, behind a discoverable constructor.

But "named solution" covers three structurally different cases, and
only one of them is cleanly shippable from the in-tree references.  A
reconnaissance of every reference paper (recorded in worklog 031)
sorted the named solutions found:

  1. **Point-IC named solutions** — the paper gives `(u, u')` at a
     specific point to high precision.  These drop straight into a
     `PainleveProblem` as the initial condition.  **Two qualify with
     gold-standard (16-digit) data:**

       - **PI tritronquée** — `u(0) ≈ -0.1875543083404949`,
         `u'(0) ≈ 0.3049055602612289`.  Source:
         `references/markdown/FW2011_painleve_methodology_JCP230/
         FW2011_painleve_methodology_JCP230.md:224-229` — computed by
         FW with an extended-precision (32-digit) Maple BVP solver
         over `[-20i, 20i]`.
       - **PII Hastings–McLeod** (`α = 0`) —
         `(u(0), u'(0)) ≈ (±0.3670615515480784, ∓0.2953721054475501)`.
         Source:
         `references/markdown/FW2014_second_PII_exploration_FoCM14/
         FW2014_second_PII_exploration_FoCM14.md:252-258`.  The two
         sign-symmetric copies "only differ in sign" (the PII `α = 0`
         symmetry `u → -u`).

  2. **Asymptotic-BC named solutions** — the paper defines the
     solution by its behaviour as `z → ±∞` (Hastings–McLeod for
     general `α`, the secondary HM, Ablowitz–Segur, the Boutroux
     `k = 0` family).  Turning an asymptotic BC into a point IC is a
     *connection problem* — a BVP solve or an asymptotic
     initialisation — not a constructor.  **Out of scope.**

  3. **Closed-form parametrised families** — exact formulas: the PII
     rational solutions `uₙ`, the PII Airy solutions for
     `α = n + ½`, the PIV entire solutions `-2z` and `-⅔z`.  These are
     genuinely shippable, but they are a *different API shape*
     (`pii_rational(n)`, `pii_airy(θ)` — parametrised families) and
     they come with *exact closed-form oracles*, which warrants its
     own self-validating-test pattern and its own ADR.  **Deferred to
     bead `padetaylor-icf`.**

So v1 is the two point-IC named transcendents.  They are also the two
*headline* ones — the tritronquée drives the whole FW 2011 figure
programme, and Hastings–McLeod is the PII solution behind the
Tracy–Widom distribution.

## Decision

Add named-transcendent **constructors** — thin functions that build
the `PainleveProblem` for a specific, literature-pinned solution and
tag it.

**File**: `src/PainleveNamed.jl`, `include`d into `module Painleve`
(the third file of that module, after `PainleveSolution.jl`; CLAUDE.md
Rule 6's ≤200-LOC budget is per-file).

**The constructors**:

```julia
tritronquee(equation::Symbol = :I; zspan = (0.0, 10.0), order = 30)
hastings_mcleod(; branch::Symbol = :positive, zspan = (0.0, 10.0), order = 30)
```

Each returns a `PainleveProblem` with (a) the cited initial condition
baked in, (b) `name` set to `:tritronquee` / `:hastings_mcleod`, and
(c) the equation parameters fixed (`hastings_mcleod` fixes `α = 0`).
The caller still chooses `zspan` (the integration window — a
placeholder for the path-network use-case) and `order`; what the
caller does *not* choose is the IC, because the IC *is* the named
solution.

`hastings_mcleod`'s `branch` selects between the two sign-symmetric
copies (`:positive` → `u(0) > 0`, `:negative` → `u(0) < 0`).

### The `name` field on `PainleveProblem`

For the named-solution tag to reach `PainleveSolution.name`, it must
travel on the `PainleveProblem` through the forwarding methods.
`PainleveProblem` therefore gains a final field
`name::Union{Symbol,Nothing}`:

  - the plain `PainleveProblem(equation; …)` constructor and every
    `_build_*` helper set `name = nothing` (a problem built from
    explicit ICs is not a named solution);
  - the named constructors set it to the solution's symbol;
  - `_painleve_solution(pp, raw)` passes `pp.name` through, so
    `solve_pade(tritronquee(:I))` yields a `PainleveSolution` whose
    `solutionname` is `:tritronquee`.

ADR-0007's `_painleve_solution` already carried a `name` parameter for
exactly this; this ADR wires the source of that argument.

## Fail-loud behaviour (CLAUDE.md Rule 1)

  - `tritronquee` with `equation ≠ :I` → `ArgumentError`: only PI has
    an in-tree tritronquée IC; the PII / PIV "tronquée" solutions in
    the references are asymptotic-BC-defined and out of scope.
  - `hastings_mcleod` with `branch ∉ (:positive, :negative)` →
    `ArgumentError` naming the two valid branches.
  - Either constructor with `zspan[1]` not equal to the IC point the
    named values are defined at (`0.0` for both) → `ArgumentError`:
    the cited `(u, u')` are the solution's value *at that point*;
    placing them at a different `zspan[1]` would silently define a
    different IVP.  The constructor will not paper over the mismatch.

## Scope

**v1 (bead `padetaylor-c86`)**:
  - `src/PainleveNamed.jl` — `tritronquee`, `hastings_mcleod`.
  - `name` field on `PainleveProblem`, plumbed into `_painleve_solution`.
  - Re-export `tritronquee`, `hastings_mcleod` from the umbrella.
  - Tests: the constructed ICs equal the cited reference values
    verbatim; `name` propagates through `solve_pade` /
    `path_network_solve` into `solutionname(sol)`; the HM `:positive`
    and `:negative` branch solutions are negatives of each other (the
    PII `α = 0` symmetry — a genuine cross-check, not a tautology);
    the fail-loud guards.  Mutation-prove (perturb a stored IC digit →
    RED; drop the `name` plumbing → RED).

**Deferred**:
  - `padetaylor-icf` — closed-form named families (PII rational /
    Airy, PIV entire) with exact-oracle self-validation.
  - Asymptotic-BC named solutions (general-`α` HM, secondary HM,
    Ablowitz–Segur, Boutroux `k = 0`) — these need a connection-problem
    solve; a constructor cannot produce them from a point IC.  If
    pursued, the right shape is a `hastings_mcleod(α)` that runs a BVP
    internally, which is a meaningfully larger piece of work.

## Rejected alternatives

  1. **A `painleve_solution(:tritronquee)` value-function returning
     numbers** — rejected for the same reason ADR-0006 rejected a
     `painleve()` value-function: the package's job is to *solve*, and
     the named constructor's job is to *pose the problem*.  It returns
     a `PainleveProblem`; the user solves it with whichever solver
     fits their domain.

  2. **Named constructors that bypass the solver and return the
     closed form** — rejected for v1's two solutions (they have *no*
     closed form — they are BVP-computed).  This *is* the right shape
     for the deferred closed-form families (bead `padetaylor-icf`),
     where it doubles as an exact test oracle; but even there the
     constructor should still return a `PainleveProblem` so the API
     stays uniform, with the closed form used only to *derive the IC*
     and to *validate* the solve.

  3. **Defaulting the IC point away from the user entirely**
     (`tritronquee(; zend = 10.0)`) — rejected: it breaks the
     `PainleveProblem` `zspan` convention every other constructor
     uses.  Keeping `zspan` with a validated `zspan[1]` is more
     consistent and still fails loud on misuse.

  4. **Storing `name` on `PainleveSolution` only, set by the
     forwarding methods via a lookup** — rejected: the forwarding
     methods are equation-agnostic; they cannot know a problem is "the
     tritronquée" unless the problem says so.  `name` belongs on the
     problem, set at construction.

## Consequences

  - **Positive**: the two headline Painlevé transcendents become
    one-liners — `solve_pade(tritronquee(:I))`,
    `path_network_solve(hastings_mcleod(), grid)` — and the resulting
    `PainleveSolution` prints and reports its own provenance
    (`solutionname(sol) === :tritronquee`).
  - **Positive**: thin and additive — one new field on
    `PainleveProblem`, one new file; the solver substrate and the
    `PainleveSolution` wrapper are untouched.
  - **Positive**: the IC values are now *cited, tested reference
    constants* in `src/`, not folklore copied between test files
    (`phase9_tritronquee_test.jl`, `fw_fig_41_test.jl`, and
    `lattice_dispatcher_test.jl` each carry the tritronquée IC as a
    literal today).
  - **Negative / honest**: v1 covers two of the dozens of named
    Painlevé solutions.  The honest reason — most are asymptotic-BC
    defined and need a connection-problem solve — is recorded above
    and in worklog 031.  This is a *foundation* for named solutions,
    not a complete catalogue.
  - **Negative / honest**: the `name` field makes `PainleveProblem` a
    7-field struct.  Acceptable — `name` is genuinely part of the
    problem's identity once a named constructor sets it, and it
    defaults to `nothing` for the common case.

## References

  - **ADR-0006** — the `PainleveProblem` builder; its "Deferred"
    section is what this ADR implements.
  - **ADR-0007** — `PainleveSolution`; its `name` field is the
    destination of the tag this ADR sources.
  - **CLAUDE.md Rule 1** — fail loud; the three constructor guards.
  - **`references/markdown/FW2011_painleve_methodology_JCP230/
    FW2011_painleve_methodology_JCP230.md:224-229`** — the PI
    tritronquée IC, verbatim.
  - **`references/markdown/FW2014_second_PII_exploration_FoCM14/
    FW2014_second_PII_exploration_FoCM14.md:252-258`** — the PII
    Hastings–McLeod IC, verbatim.
  - **Worklog 031** — the in-tree reference reconnaissance that
    sorted the named solutions into the three cases above.
  - **Bead `padetaylor-icf`** — the deferred closed-form families.
