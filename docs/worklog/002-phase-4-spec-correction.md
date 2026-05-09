# Worklog 002 — Phase 4 (StepControl) + DESIGN/HANDOFF spec correction

**Date**: 2026-05-09
**Author**: Claude (orchestrating); opus subagent (implementation)
**Scope**: Phase 4 of Stage 2 implementation plan — `StepControl` module
plus a paper-grounded correction to DESIGN.md §4 and HANDOFF.md's
Phase-4 section.

> The take-home from this shard: **trust the paper line range, not the
> spec doc**. The spec was wrong because the previous author wrote it
> from memory of "Jorba-Zou + safety factor", not from reading the
> paper.  Triangulating against multiple ground truths (paper, TI.jl,
> mpmath, wolframscript) catches this at zero cost beyond the few
> minutes to set up the four scripts.

## What changed

Two new public functions in `src/StepControl.jl` (~111 LOC body + 128-
line literate top-of-file docstring):

- `step_jorba_zou(coefs, eps_abs; eps_rel = eps_abs)` — Jorba-Zou
  2005 §3.3.1 first step-size control, fixed-order form. Ported
  verbatim from `external/TaylorIntegration.jl/src/integrator/
  stepsize.jl:17-89`, with TI.jl's ε-resolution and the `_second_
  stepsize` fallback for all-zero candidates.
- `step_pade_root(P, z_current, target)` — FW 2011 §3.1 forward-
  projection-onto-direction heuristic. Computes denominator roots via
  `Polynomials.roots`, projects each onto the unit step direction,
  filters forward (positive projection), caps at `|target -
  z_current|`.

4 testsets / 7 individual @test calls in `test/stepcontrol_test.jl`.
Pinned oracle in `test/_oracle_stepcontrol.jl`. Total suite: 186 prior +
6 new pass = 192/192 GREEN.

## The DESIGN.md / HANDOFF.md spec correction

The original DESIGN.md §4 Phase 4 (and HANDOFF.md's mirror) specified
the Jorba-Zou step formula as

```
ρ = min((eps / |c[p]|)^(1/p), (eps / |c[p-1]|)^(1/(p-1)))
h = (ρ / e²) · exp(-0.7 / (p - 1))
```

with attribution "Jorba-Zou 2005 §3.2 eq. 3-8". This is wrong on three
counts:

1. **The `0.7` constant has no source.** Grep `0.7` and `safety` in
   both `references/markdown/JorbaZou2005_*/JorbaZou2005_*.md` (the
   paper) and `external/TaylorIntegration.jl/src/`: zero matches in
   either. The constant was a hallucination at spec-write time. No
   other Taylor-IVP impl I could find (chebfun, TI.jl) has a
   comparable safety factor. RESEARCH.md §3.5 also cites "Jorba-Zou
   §3.2.1 eq. 3-8" but that section / equation does not exist in the
   paper.

2. **The `/e²` is misplaced for fixed-`p` use.** The paper's actual
   formula (§3.3.1 eq. 11) reads
       ρ̄_j = |x[j]|^(-1/j),
       h_n = ρ_n / e²,
   where the `/e²` only makes sense at the *adaptive* `p_n =
   -½ ln ε + 1` from eq. 11. PadeTaylor fixes `p = 30` (FW 2011), so
   the `/e²` is implicitly absorbed into the `ε`-substitution that
   eliminates `ρ̄`'s dimensional dependence on the analyticity radius.
   At fixed `p`, the equivalent formula is TI.jl's
       h = min over k ∈ {p-1, p} of (ε / |c[k+1]|)^(1/k)
   (no `/e²`, no safety factor).

3. **The two indices `(p, p-1)` were correct, but the order of `min`
   is over `k`, not over `(eps/|c[p]|)^(1/p)` vs `(eps/|c[p-1]|)^(1/(p-1))`
   as separate quantities.** The mathematical intent of DESIGN.md's
   formula was right; only the safety-factor and `/e²` were spurious.

### Triangulation methodology

To prove the spec correction wasn't itself a fresh hallucination, I
captured the canonical test case
(`coefs = [1/k! for k=0..30]`, `ε = 1e-12`) in three independent ways
and checked they all give the same number:

- **Julia**: `TaylorIntegration.stepsize(Taylor1(coefs, 30), 1e-12)`
  → `4.501206370338986`
- **Python**: `mpmath` at 50 dps, formula
  `min((eps · k!)^(1/k) for k in {29, 30})`
  → `4.50120637033898607690318848315848021108523516609`
- **Mathematica**: `wolframscript` at 50 dps, same formula
  → `4.50120637033898607690318848315848021108523516609`

47 decimal digits byte-identical between mpmath and wolframscript;
TI.jl matches to all available Float64 digits. The DESIGN.md formula
gives `0.633` for the same input — clearly different. Three-source
consensus says: TI.jl-exact is right, DESIGN.md was wrong.

The capture/verification scripts live at
`external/probes/stepcontrol-oracle/{capture.jl, capture.py,
capture.wl, verify.jl}`. Re-running them is the single command
`julia --project=external/probes/stepcontrol-oracle external/probes/
stepcontrol-oracle/verify.jl` (after wolframscript + python3 captures).

### Documents updated

- `DESIGN.md §4 Phase 4` — formula corrected with explicit
  "CORRECTED 2026-05-09" note and the 47-digit consensus value.
- `HANDOFF.md` Phase-4 section — replaced the stale "Goal / Ground
  truth / Test plan / Key gotcha" subsections with a "Phase 4
  SHIPPED" status block pointing at the corrected formula.
- `src/StepControl.jl`'s top docstring — full triangulation argument
  for any future agent reading the module.

## Mutation-proof procedure (verified before commit)

Both load-bearing functions independently mutation-tested:

- **Mutation A** — off-by-one in the Jorba-Zou exponent
  (`(ε/|c[k+1]|)^(1/k) → (ε/|c[k+1]|)^(1/(k+1))`):
  fails 4.1.1 + 4.1.2 (190 pass / 2 fail).
- **Mutation C** — projection axis swap in `step_pade_root`
  (`real(...)` → `imag(...)`):
  fails 4.1.3 + 4.1.4 (190 pass / 2 fail). The real-axis paths in
  both tests give imag projection = 0, so the function returns the
  uncapped distance instead of the correct 2.0 / 3.0.

Both mutations were applied, observed to RED, then restored before
the commit.

## Frictions surfaced (for the next agent)

### 1. The TaylorIntegration.jl UUID has changed.

When setting up `external/probes/stepcontrol-oracle/Project.toml`, I
copied the UUID from somewhere (probably an old DESIGN.md draft) and
got `92b13dbe-c966-51a2-8445-caca9f8c7bd4`. The current registered UUID
is `92b13dbe-c966-51a2-8445-caca9f8a7d42` (different last 12 hex
digits). `Pkg.add` errors with "expected package … to be registered;
you may have provided the wrong UUID". When adding TI.jl to a probe
env, **either omit the UUID and let Pkg resolve, or copy the UUID from
the registered General entry**.

### 2. `Polynomials.roots` returns ComplexF64 even for real input.

`Polynomials.roots(Polynomial([1.0, -0.5]))` returns
`ComplexF64[2.0 + 0.0im]`, not `[2.0]`. The downstream projection
arithmetic must handle complex arithmetic naturally; the `real(...)`
extraction picks up the right value. This worked correctly in the
shipped impl but caught me by surprise during oracle generation.

### 3. Arb-Polynomial root-finding is still open.

`Polynomials.roots ∘ Polynomials.Polynomial` is unvalidated for
`T = Arblib.Arb`. Per HANDOFF.md, this was already flagged as a
v2-deferred item.  I did not probe it in Phase 4 because the Phase-4
test corpus is `Float64`-only; `step_pade_root` for Arb Padé
approximants will need a separate companion-matrix root finder via
`GenericSchur.jl` or the Arb library's own root finder. File a bead
when Phase 8 (`PadeTaylorArblibExt`) actually exercises this path.

### 4. Mid-script `cd` persists across Bash tool calls.

I `cd`'d into `external/probes/stepcontrol-oracle/` for one capture
script run, then in the next tool call expected the cwd to be the
project root.  It wasn't — `cd` persists between calls.  **Use
absolute paths in long-running tool sessions** (or always re-`cd`
explicitly).

## Acceptance

| Phase | Module | LOC | Tests | Mutation-proof | Acceptance |
|---|---|---|---|---|---|
| 4 | `StepControl` | 111 | 4 testsets / 7 @test | ✅ A bites 4.1.1+4.1.2; C bites 4.1.3+4.1.4 | 192/192 GREEN |

## Pointers

- [`src/StepControl.jl`](../../src/StepControl.jl) — the module + the
  triangulation argument in its top docstring.
- [`test/stepcontrol_test.jl`](../../test/stepcontrol_test.jl) — the
  4 testsets + verified mutation-proof block.
- [`external/probes/stepcontrol-oracle/`](../../external/probes/stepcontrol-oracle/)
  — capture scripts (Julia / Python / wolframscript) + verifier.
- [`docs/worklog/001-stages-Z-1-2-handoff.md`](001-stages-Z-1-2-handoff.md)
  — the previous shard, for context.
