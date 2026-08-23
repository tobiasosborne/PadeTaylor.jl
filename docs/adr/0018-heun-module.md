# ADR-0018 — Heun module: HeunG + HeunC as first-class citizens

**Status**: Accepted (2026-05-18) | **Bead**: epic `padetaylor-72w` (Heun
arc) | **Worklog**: 051 (forthcoming)

## Decision

Ship a new `src/Heun.jl` module that evaluates the general Heun function
`HeunG(a, q, α, β, γ, δ; z)` and the confluent Heun function
`HeunC(q, α, γ, δ, ε; z)` at arbitrary complex `z`, by:

1. Computing the Frobenius power series at `z = 0` via the
   DLMF 31.3.2–31.3.4 three-term recurrence (HeunG) and the
   Motygin 2018 eqs. 3–4 three-term recurrence (HeunC);
2. Sampling the series at a small `z₀` (default `ε = 0.05`) to obtain
   `(u(z₀), u'(z₀))` as initial conditions for `PadeTaylorProblem`;
3. Walking the path-network from `z₀` to the user's target `z`,
   reusing all existing path-network + branch-tracker machinery
   (FW 2011 §3.1 wedge tree, BranchTracker for the regular singular
   points at `z = 0, 1, a`).

Public API in `src/Heun.jl`:

```julia
heun_g(a, q, α, β, γ, δ; z, h=0.1, order=30, sheet=:principal, ...) -> Complex
heun_c(q, α, γ, δ, ε; z, h=0.1, order=30, sheet=:principal, ...) -> Complex
```

with element-type promotion (`Float64` or `BigFloat`) and the full
`path_network_solve` kwarg surface forwarded through.  Both functions
return `1.0 + 0im` at `z = 0` by normalisation.

Adopt **DLMF §31.2 / §31.12 parameter conventions** as the canonical
form (matches Mathematica `HeunG[a,q,α,β,γ,δ,z]` and `HeunC[q,α,γ,δ,ε,z]`
exactly; matches Maple HeunG; diverges from Maple HeunC where we document
the mismatch).

Branch handling: default `cross_branch = false` with
`branch_points = (0, 1, a)` for HeunG and `(0, 1)` for HeunC, and
`branch_cut_angles = π` per branch (Wolfram convention).  Opt-in
sheet bookkeeping via `sheet` kwarg.

## Context

Heun's general equation (DLMF 31.2.1) is the next layer up from the
hypergeometric ODE: instead of three regular singular points it has
four (at `z = 0, 1, a, ∞`), with Fuchsian constraint
`α + β + 1 = γ + δ + ε`.  Confluent Heun (DLMF 31.12.1) merges two
of those into an irregular singularity at infinity.

These functions appear in: black hole quasinormal modes (Teukolsky
master equation reduces to HeunC, Fiziev 2009
`references/markdown/heun/teukolsky/`), the hydrogen molecule ion
H₂⁺, Bethe ansatz eigenfunctions in integrable spin chains, two-centre
Coulomb problems, dislocation theory.  No good open-source numerical
implementation has existed.

Existing tools:

  - **Mathematica `HeunG` / `HeunC`** (built-in since v9, ~2012):
    closed-source; documented as fragile past the second finite
    singularity; ground-truth quality only inside the unit disc and
    on small extensions.  Wolfram blog 2020-05-06: *"We do not know
    their explicit forms and are forced to work with their generating
    equations."*
  - **Maple `HeunG` / `HeunC`** (built-in): Motygin 2015 abstract:
    *"the only software package able to evaluate [general] Heun
    functions numerically is Maple™… [its] implementation is imperfect"*
    and *"code is not available."*
  - **Motygin 2015 / 2018 / 2020** (`external/heun-refs/motygin-Heun/`
    and `motygin-confluent-Heun/`): GPL-3.0 MATLAB/Octave packages
    using overlapping-circle analytic continuation.  Actively maintained
    (last update 2025-10-27 for confluent).  **The current state of
    the art in open-source Heun numerics.**

PadeTaylor.jl already implements every primitive we need:

  - Local Taylor → Padé conversion (Phase 2, `RobustPade`);
  - Path-network walking with wedge candidates around singularities
    (`PathNetwork.path_network_solve`);
  - BranchTracker for sheet bookkeeping across cuts (ADR-0013, worklog
    042);
  - Adaptive FFW step controller (ADR-0011);
  - Frobenius-series IC pattern (already used for Painlevé asymptotic ICs).

The Heun module is therefore not a new algorithm; it is the *systematic
application* of the path-network machinery to a different class of ODE.

## Architecture

### Module layout

```
src/Heun.jl                    -- HeunG + HeunC public API + Frobenius IC
                                  factories.  ~250–300 LOC target; split
                                  to src/HeunC.jl if it grows past Rule 6's
                                  200 effective LOC.
src/PadeTaylor.jl              -- include and re-export heun_g, heun_c
test/heun_test.jl              -- golden-master suite (Stage 5)
external/probes/heun-oracle/   -- wolframscript oracle (Stage 2 — shipped)
external/heun-refs/            -- Motygin reference impl (Stage 1 — shipped)
docs/heun_research.md          -- background synthesis (Stage 0 — shipped)
docs/heun_dlmf_summary.md      -- recurrence reference (Stage 0b — shipped)
docs/teukolsky_heun_mapping.md -- Teukolsky → HeunC parameter map
                                  (Stage 2d — shipped)
docs/adr/0018-heun-module.md   -- this file
docs/worklog/051-heun-functions.md -- session log (Stage 6)
```

### The HeunG Frobenius recurrence

DLMF 31.3.2–31.3.4 (verified against `references/heun/dlmf-31/31.3.html`
and against the wolframscript oracle at
`external/probes/heun-oracle/oracles.txt`):

```
c_0 = 1
c_1 = q / (a γ)

For j ≥ 1:
  R_j · c_{j+1} = (Q_j + q) · c_j - P_j · c_{j-1}

  P_j = (j - 1 + α)(j - 1 + β)
  Q_j = j · ((j - 1 + γ)(1 + a) + a δ + ε)
  R_j = a (j + 1) (j + γ)

  ε = α + β + 1 - γ - δ          (Fuchsian)
```

Radius of convergence: `min(1, |a|)` (distance to nearest finite
singularity from origin).

### The HeunC Frobenius recurrence

Motygin 2018 eqs. 3–4 (cited in `docs/heun_research.md:563-579`):

```
b_0 = 1
b_{-1} = 0

For n ≥ 1:
  P_n · b_n = Q_n · b_{n-1} + R_n · b_{n-2}

  P_n = n (n + γ - 1)
  Q_n = -q + (n - 1)(γ + δ - ε + n - 2)
  R_n = α + (n - 2) ε
```

Radius of convergence: 1 (nearest singular point z = 1).  Initial values:
`b_0 = 1`, and `b_1 = Q_1·b_0 / P_1 = (-q + 0) / γ = -q/γ`.

(Note: DLMF Ch. 31 does NOT publish the confluent-Heun recurrence
explicitly; it refers the reader to Ronveaux 1995 Part B.  We adopt
Motygin's form — algebraically the cleanest available — and validate
it against the wolframscript oracle.)

### Path-network IC delivery

To use `path_network_solve` we need an IC tuple `(u(z₀), u'(z₀))` at a
finite starting point `z₀ ≠ 0`.  Strategy:

1. Choose `z₀ ∈ (0, ε_start]` real with `ε_start = 0.05` (default;
   well inside the Frobenius disc and far enough from `z = 0`'s
   coefficient pole to avoid conditioning issues).
2. Evaluate the Frobenius series at `z₀` truncated at `N_taylor = 60`
   terms (FFW 2017 §2.1.1 default).
3. Differentiate term-by-term to evaluate `u'(z₀)`.
4. Pass `(u(z₀), u'(z₀))` as `y0` to `PadeTaylorProblem`.
5. For the path-network walker, set `zspan = (z₀, z_target)` and
   `branch_points = (0+0im, 1+0im, complex(a))` for HeunG (omit `a`
   for HeunC).

When the user requests `z` with `|z| < ε_start`, return the Frobenius
series sum directly — no path-network walk needed.

### RHS factories

Both Heun equations can be written as `u'' = f(z, u, u')`:

**HeunG** (DLMF 31.2.1, rearranged):
```
u'' = -((γ/z + δ/(z-1) + ε/(z-a)) · u'
        + (α β z - q) / (z(z-1)(z-a)) · u)
```

**HeunC** (DLMF 31.12.1, rearranged):
```
u'' = -((γ/z + δ/(z-1) + ε) · u'
        + (α z - q) / (z(z-1)) · u)
```

These are thread-safe closures (the path-network walker runs the five
wedge candidates in parallel via `Threads.@threads` in
`_wedge_evaluations` — `f` must capture no mutable state).

### Branch handling at the regular singular points

`z = 0, 1, a` for HeunG (and `z = 0, 1` for HeunC) are *fixed*
singularities of the linear ODE — not movable poles like Painlevé.
The path-network walker treats them by *avoidance*: the BranchTracker
refuse-mode (default `cross_branch = false`) rejects any wedge
candidate whose step segment would cross a branch cut.  Concretely
this means the walker detours into the upper or lower half-plane
when the goal lies on the other side of `z = 1` or `z = a`.

For the cross-branch mode (Riemann-surface monodromy queries),
`cross_branch = true` is exposed as a kwarg, but with v1 we don't
attempt to validate the resulting sheet bookkeeping against
Mathematica's monodromy choices — this is a known gap, deferred to a
follow-up bead (probably `padetaylor-?heun-monodromy`).

### Why not a separate Frobenius recurrence for HeunC

DLMF 31.12 doesn't publish an explicit recurrence.  Motygin 2018 does.
We adopt Motygin's form because:

1. It's the published recurrence used by the only actively-maintained
   open-source HeunC implementation.
2. It's algebraically clean (3 terms per coefficient, no special cases
   except `b_{-1} = 0`).
3. It validates against Mathematica's HeunC values to ≥ 10 digits on
   the cross-validation oracle (regime A of
   `external/probes/heun-oracle/oracles.txt`).

We *could* instead skip the Frobenius series entirely and integrate
DLMF 31.12.1 as an ODE from `z = ε_start` with `u(0) = 1`, `u'(0) =
-q/γ` extrapolated.  But the Frobenius approach at `z₀` is
(a) faster (a recurrence is O(N), the ODE walk is O(steps × order²)),
(b) more accurate when `|z₀| < ε_start ≪ 1`, and (c) gives us a
free unit test (the recurrence values must match path-network values
in the overlap region).

## Locked decisions

1. **Two functions, one module.** `heun_g` and `heun_c` ship together
   in `src/Heun.jl`.  Split only if effective LOC exceeds 200
   (CLAUDE.md Rule 6).
2. **DLMF parameter conventions.** `heun_g(a, q, α, β, γ, δ; z)` and
   `heun_c(q, α, γ, δ, ε; z)`.  No Maple-style aliases at v1; document
   the convention difference in module docstring.
3. **Frobenius-series IC at `z₀ = 0.05`.** Configurable via kwarg
   `epsilon_start`.  N=60 Taylor terms (FFW 2017 default).
4. **Element type follows the user.** Same promotion rules as
   `PadeTaylorProblem`; `BigFloat`-256 fallback available by passing
   `BigFloat`-typed parameters.
5. **Branch handling: refuse mode default** (`cross_branch = false`).
   Cross-branch is exposed but its Mathematica-monodromy-equivalence
   is unvalidated at v1.
6. **Oracle pinning.** Tests pin against the wolframscript oracle at
   `external/probes/heun-oracle/oracles.txt`.  Two-oracle validation
   (Mathematica + Motygin MATLAB) is in flight (Stage 2c, bead
   `padetaylor-3it` cross-val task).
7. **Wall-time gate: 2 s per test** (user mandate 2026-05-18 in this
   conversation).  Each HG.* / HC.* / TK.* test must complete in
   under 2 s; full `test/heun_test.jl` under 60 s.

## Alternatives considered

1. **Hypergeometric reduction**.  When `q` takes special "apparent
   singularity" values, HeunG degenerates to ₂F₁ products (DLMF
   §31.7).  We could special-case these.  *Rejected for v1:* the
   apparent-singularity values are themselves discrete; detecting
   them numerically is error-prone; and the path-network approach
   handles them as well as any other parameter set.

2. **Pure Frobenius + analytic continuation** (à la Motygin 2015).
   Use overlapping circles with explicit re-expansion at each circle
   centre.  *Rejected* because we already have the path-network
   machinery; using it gives us automatic adaptive step control,
   branch tracking, and Stage-2 dense interpolation for free.

3. **Mathematica wrapper**.  Shell out to wolframscript for every
   evaluation.  *Rejected*: defeats the purpose of having a Julia
   implementation, has 100 ms+ startup overhead per call,
   non-portable.

4. **Use AAA rational approximation** (Nakatsukasa-Sète-Trefethen
   2018) once we have sample values.  *Out of scope for v1* but an
   interesting future direction: build AAA on a dense path-network
   sample set to get a global rational approximation usable everywhere
   without re-walking.

## Consequences

  - PadeTaylor.jl becomes (to our knowledge) the **first open-source
    Julia implementation of Heun functions** — and the first
    open-source non-MATLAB implementation period (Motygin is
    GPL-3.0 MATLAB).
  - The Heun module exercises BranchTracker on its first non-Painlevé
    target, validating that A4's design (worklog 042) was sound.
  - Wall-time discipline (2 s per test, user mandate) forces early
    perf attention.  If single-point HeunG evaluation comes in at
    > 1 s, we'll need either lower-order Padé (15-15 vs 30-30) or a
    caching layer for repeat-parameter calls.

## Open follow-ups (deferred beads)

  - **Heun monodromy validation**: confirm `cross_branch = true` sheet
    bookkeeping matches Mathematica's `MeijerG`-style monodromy choice.
  - **Logarithmic-γ regularisation**: Motygin 2020's `Hl☼` function
    that's C∞ in `γ` across integer values.  Useful for QNM problems
    where `γ` crosses integer boundaries during parameter sweeps.
  - **HeunB / HeunD / HeunT** (biconfluent / doubly-confluent /
    triconfluent).  Specialised physics applications.
  - **AAA-global wrapper**: after the path-network builds a dense
    sample set, fit a global AAA rational for repeated evaluations.

## References

  - **DLMF Chapter 31** — `references/heun/dlmf-31/{31.2,31.3,31.4,31.12}.html`.
  - **`docs/heun_dlmf_summary.md`** — line-cited DLMF extracts; the
    authoritative recurrence transcription.
  - **`docs/heun_research.md`** — broader survey including Motygin's
    HeunC recurrence (line 563–579).
  - **`docs/teukolsky_heun_mapping.md`** — Fiziev 2009 Eq. II.6
    parameter map for the Stage 4c QNM demo.
  - **`external/probes/heun-oracle/oracles.txt`** — 42-record
    wolframscript oracle, the test-suite ground truth.
  - **`external/heun-refs/README.md`** — implementation inventory;
    Motygin's MATLAB code is the cross-validation oracle.
  - **ADR-0004** — path-network architecture.
  - **ADR-0011** — adaptive FFW step controller.
  - **ADR-0013** — constrained-wedge routing + BranchTracker.
  - **CLAUDE.md Laws 1 & 2, Rules 6 & 10.**
