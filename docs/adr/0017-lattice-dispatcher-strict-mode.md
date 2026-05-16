# ADR-0017 — Lattice-dispatcher: edge-gated IVP source + strict/fail-soft mode

**Status**: Accepted (2026-05-16) | **Bead**: `padetaylor-0tj` | **Worklog**: 032 (discovery).

## Decision

`LatticeDispatcher.lattice_dispatch_solve` v3 makes two paired
changes: (a) when no `mask` is supplied, Step 1 routes the IVP
through `EdgeGatedSolve.edge_gated_pole_field_solve` instead of plain
`PathNetwork.path_network_solve`, so the BCs the per-row BVP reads
come from genuine pole-field cells rather than IVP-corrupted smooth
cells; and (b) a new `strict::Bool = true` kwarg gates the response
to per-row BVP non-convergence — `strict = true` preserves the v1/v2
fail-fast throw, `strict = false` catches the exact "Newton did not
converge" `ErrorException`, tags the affected cells `:bvp_fail`, and
continues.  Default `strict = true` keeps every existing test
invariant; the edge-gated IVP source becomes the new default because
the old path was *wrong* (not merely slow) for the documented
PI tritronquée use case.

## Context

Bead `padetaylor-0tj`: on a 121×121 lattice over `[-20,20]²`, PI
tritronquée's pole-free angular sectors return FW2011_*.md:401
-corrupted values —

> "smooth regions are unstable regions for any IVP solver, and the
> associated loss of accuracy ought not be carried back into the pole
> field. One way to prevent this is to force the path selection
> algorithm to complete pole fields before stepping into smooth
> regions … this process can readily be fully automated."
> — `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:401`

Step 3 then read those corrupted cells as BVP Dirichlet BCs, Newton
diverged on the first non-convergent row, and `bvp_solve` threw
`ErrorException("bvp_solve: Newton did not converge…")`.  The throw
was Rule-1-correct fail-loud; the upstream input was the bug.
Worklog 032 traced the corruption to the plain `path_network_solve`
at the old `src/LatticeDispatcher.jl:206`.  `EdgeGatedSolve` is FW
md:401 made executable and is the right cure.

## Architecture

### (a) Edge-gated IVP source (default, `mask === nothing`)

Step 1 calls `edge_gated_pole_field_solve(prob, xs, ys; h = h_path,
order, edge_level)`, takes `u_grid = gated.u_grid` (NaN outside the
gated field), rebuilds `up_grid` by scattering
`gated.pn_solution.grid_up` onto the one-ring-dilated-field cells,
and sets the partition `mask_used = gated.field_mask` (the
un-dilated edge-confirmed mask — only edge-confirmed pole-field
cells serve as BVP flanks; one-ring frontier cells are smooth-side).

When `mask !== nothing` (FW md:401 manual workflow, documented at
`LatticeDispatcher.jl:81-83` since v1), Step 1 falls back to plain
`path_network_solve` over the full grid and the caller's mask
drives Step 3 — preserving the manual-classification API contract.

### (b) `strict::Bool = true` kwarg

Narrowly-scoped `try/catch` around the `bvp_solve` call:

```julia
catch err
    if !strict && err isa ErrorException &&
       occursin("bvp_solve: Newton did not converge", err.msg)
        for k in run_start:run_end; region_tag[k, j] = :bvp_fail; end
        nothing
    else
        rethrow()
    end
end
```

Matches on the exact message from `src/BVP.jl:332-338`.  Any other
exception (`ArgumentError`, `InexactError`, etc.) rethrows.

`:bvp_fail` joins `{:ivp, :bvp, :ivp_only}` in
`LatticeSolution.region_tag`'s enum.

## Alternatives considered, rejected

- **Keep plain `path_network_solve`, just add `strict`.**  Rejected:
  BCs still corrupted; `strict = false` would mask wrong answers
  (plausible BVP fills from corrupted BCs) — Rule 1 violation by
  stealth.
- **Always edge-gated, no manual fallback.**  Rejected: `mask` kwarg
  is documented v1 public API (`LatticeDispatcher.jl:81-83`); a caller
  who supplied a classification has committed to it and silent
  override would be a semantic break.
- **Auto-retry with larger `N_bvp` on non-convergence.**  Rejected
  for v1: orthogonal to the BC-corruption fix and not actionable
  without knowing whether non-convergence is "needs resolution"
  (helpable) vs "no solution exists" (not).  Follow-up bead.
- **One-ring-dilate `field_mask` as the partition.**  Rejected:
  frontier cells weren't edge-confirmed, just solved opportunistically.
  Using them as BVP flanks re-introduces the kind of stealth wrongness
  this ADR removes.

## Consequences

- Smooth cells outside the gated field are `NaN + NaN·im` by default
  — Rule 1 honest; old code relying on a value at every cell now
  sees NaN unless the per-row BVP fills it.
- `strict = true` default preserves every existing test invariant
  byte-for-byte (back-compat).
- `:bvp_fail` is a new region tag; downstream consumers pattern-
  matching on the enum need a fifth arm (none currently do).
- Manual-classification (`mask=`) workflow unchanged.

### Test coverage (separate sub-agent task)

Target `test/lattice_strict_test.jl`:
  - Default `strict = true` still throws (back-compat guard).
  - `strict = false` yields `:bvp_fail` and no throw on the same
    fixture.
  - PI tritronquée on coarse `[-20,20]²` completes under
    `strict = false`; smooth-sector cells are NaN.
  - `mask !== nothing` fallback still uses plain `path_network_solve`.
  - Mutation-proof: perturb the `occursin` substring → still passes
    strict-throw (catch too narrow); perturb to one-ring partition →
    flank misclassification detected.

### Open follow-ups

- Auto-retry-with-larger-N_bvp on non-convergence.
- A `:bvp_fail` row's residual / message is currently discarded; a
  diagnostics vector on `LatticeSolution` could surface it.

## References

- Bead `padetaylor-0tj` — symptom and locked fix spec.
- FW 2011 md:401 — motivation; the "smooth regions are unstable" line.
- Worklog 032 — discovery of the BC-corruption failure mode.
- ADR-0004 — `docs/adr/0004-path-network-architecture.md`.
- `src/EdgeGatedSolve.jl` — the gated solver adopted as default
  IVP source.
- `src/BVP.jl:332-338` — exact `ErrorException` whose message the
  fail-soft catch matches.
- CLAUDE.md Rules 1, 6, 10.
