# Worklog 026 — `PainleveProblem` layer + `prob.order` footgun fix

**Date**: 2026-05-14 (continues worklogs 022–025)
**Author**: Claude Opus
**Scope**: Two beads, done as one arc: `padetaylor-9xf` (the
`path_network_solve` / dispatcher `prob.order` footgun) first, because
ADR-0006 flagged that the new layer would otherwise have to work
around it; then `padetaylor-avl` — implement the `PainleveProblem`
builder layer per ADR-0006.

> **Take-home**: `src/Painleve.jl` ships the per-equation problem
> builder — `PainleveProblem(:I…:VI; <params>, u0, up0, zspan)` —
> assembling the correct `PadeTaylorProblem` for any of the six
> Painlevé equations: RHS selection (3 new direct-frame factories for
> PI/PII/PIV; PIII/PV/PVI wrap the existing `CoordTransforms` /
> `SheetTracker` factories), coordinate-frame handling, parameter
> bookkeeping, branch-point guards. It is a *builder*, not a
> `painleve()` value-function (ADR-0006's load-bearing decision). Test
> suite **1432 → 1487 GREEN** (+7 from the 9xf regression tests, +48
> from `painleve_test.jl`'s 47 assertions plus the umbrella-load
> check). All mutation-proofs bite.

## Part 1 — `padetaylor-9xf`: drivers now honour `prob.order`

`path_network_solve`, `dispatch_solve`, and `lattice_dispatch_solve`
each carried their own `order::Integer = 30` kwarg and silently
ignored the `order` field of the `PadeTaylorProblem` they were handed
— so a problem built at any order ≠ 30 was solved at 30 regardless.
Caught while writing `figures/fw2011_fig_5_2.jl` (worklog 025): an
`(order, h)` sweep that was accidentally all-order-30.

**Fix**: default the kwarg to `prob.order` (resp. `prob_ivp.order`).
An explicit `order =` still overrides. Behaviour-preserving for the
whole pre-existing suite — every test either builds its problem at
order 30 (new default == old default) or passes `order` explicitly.
Note that fixing `path_network_solve` alone already fixes
`dispatch_solve`'s internal path-network calls by propagation
(dispatch builds its segment problems with `order`, then hands them
down); `dispatch_solve`'s own kwarg default was changed too for
consistency.

**Regression tests**, all mutation-proven (revert `= prob.order` →
`= 30`, confirm RED, restore):
  - `PN.7.1` (`pathnetwork_test.jl`) — default-order solve == explicit
    `prob.order` solve, ≠ explicit order-30 solve.
  - `DP.1.3` (`dispatcher_test.jl`) — same for `dispatch_solve`.
  - `LD.2.2` (`lattice_dispatcher_test.jl`) — same for
    `lattice_dispatch_solve`.

Committed as `a2f319c`.

## Part 2 — `padetaylor-avl`: the `PainleveProblem` layer

### `src/Painleve.jl` (new)

Per ADR-0006. The module:

  - **3 new direct-frame RHS factories** — `_pI_rhs()` (`6u²+z`),
    `_pII_rhs(α)` (`2u³+zu+α`, FW 2014 eq. 2), `_pIV_rhs(α,β)`
    (Reeger–Fornberg 2014). The PIV closure has a `iszero(u)` guard
    that throws `DomainError` rather than letting an `Inf` from the
    `β/u` / `(u')²/(2u)` terms leak into the Taylor jet (PIV solutions
    carry movable zeros — CLAUDE.md Rule 1).
  - **`PainleveProblem{P,TF,FF}`** — wrapper struct carrying
    `equation`, `params`, the underlying `problem`, `frame`
    (`:direct` | `:transformed`), and the `to_frame` / `from_frame`
    coordinate maps.
  - **6 constructors** dispatched by equation symbol. PI/PII/PIV build
    a direct-frame problem; PIII/PV/PVI delegate to the *existing*
    `pIII_transformed_rhs` / `pV_transformed_rhs` /
    `pVI_transformed_rhs` factories, map the IC into the ζ-frame, and
    store the transform pair. Branch-point guards: `zspan[1]` may not
    be `z = 0` (III/V) or `z ∈ {0,1}` (VI).
  - **`_validate_kw`** — fail-loud on unknown equation, missing
    required parameter, or a parameter the equation does not take.
  - **Forwarding methods** — `solve_pade(pp; …)` and
    `path_network_solve(pp, grid; …)`. `:direct` problems forward
    transparently to `pp.problem` (clean now that 9xf is fixed);
    `:transformed` problems throw (see Frictions).

Wired into the umbrella module (`include` + `using` + `export
PainleveProblem`); `runtests.jl` gains the umbrella-load assertion and
`painleve_test.jl`.

### `test/painleve_test.jl` (new) — 47 assertions

`PV.1.*` constructor correctness (RHS sample-point match, frame,
params, ζ-mapped IC for the transformed equations); `PV.2.1`
fail-loud guards; `PV.3.1` forwarding (`:direct` ≡ direct call,
`:transformed` throws); `PV.4.1` end-to-end PII pole-field smoke test.
Mutation-proven — three mutations (perturb the PII factory, perturb
the PIV factory, drop the `_forward_guard` call) each bite their
target testset; full procedure in the file footer.

## Acceptance

ADR-0006's v1 scope: `src/Painleve.jl` + 6 constructors + 3 new RHS
factories + forwarding methods + tests + mutation-proof — all shipped.
The three implementation refinements (no separate `z0` kwarg,
`:transformed` forwarding throws rather than auto-maps, no `bvp_solve`
forwarding method) are recorded in ADR-0006's "Implementation
refinements" section. ADR-0006 status → Accepted.

## Frictions

  - **A ζ-frame solution cannot be faithfully round-tripped.** The ADR
    sketch had the forwarding methods map `:transformed` results back
    to the z-frame. But a `PathNetworkSolution`'s per-node Padé store
    and step length `h` are intrinsically ζ-frame — a Padé in ζ is not
    a Padé in z, and `h` in ζ is not constant in z. Mapping the
    *point-value* fields (`grid_*`, `visited_z/u/up`) back is clean;
    mapping the Padé store is not. Rather than return a half-mapped
    result (silent confusion, against Rule 1), v1's forwarding throws
    for `:transformed` and points the caller at `pp.problem` +
    `pp.from_frame`. Faithful point-value remapping is bead
    `padetaylor-soi` (P3).
  - **`bvp_solve` doesn't take a problem.** The ADR listed a
    `bvp_solve` forwarding method, but `bvp_solve` takes RHS closures
    + boundary data, not a `PadeTaylorProblem` — there is no
    well-typed `bvp_solve(pp)`. Dropped from v1.
  - **`path_network_solve`'s `order` kwarg was the rake ADR-0006 saw
    coming.** Part 1 cleared it before Part 2, exactly as the ADR's
    "Constraints carried from existing code" section asked — so the
    `:direct` forwarding methods are a clean one-line unwrap with no
    `order`-threading workaround.

## Beads

  - `padetaylor-9xf` — **closed** (Part 1).
  - `padetaylor-avl` — **closed** (Part 2).
  - `padetaylor-soi` **P3** (filed) — `:transformed` forwarding with
    point-value remapping.

## Hard-won lesson

**An ADR written one turn can clear its own blocker the next.**
ADR-0006's "Constraints carried from existing code" section named the
`prob.order` footgun and said "fix it first or the forwarding methods
work around it." Doing exactly that — 9xf before avl — meant the
PainleveProblem forwarding methods are a clean unwrap instead of
carrying an `order`-threading wart. Writing the constraint down in the
ADR is what made the ordering obvious.
