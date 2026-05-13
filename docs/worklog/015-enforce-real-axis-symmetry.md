# Worklog 015 — PathNetwork `enforce_real_axis_symmetry` kwarg (bead `padetaylor-dtj`)

**Date**: 2026-05-13
**Author**: Claude Opus
**Scope**: Promote the user-space Schwarz-reflection workaround from
`examples/tritronquee_3d.jl` (worklog 014 §"Bug 1") to a library-level
opt-in kwarg.  4 new testsets, 13 new assertions, mutation-proven on
both the conj mirror and the upper-canon filter.  GREEN at 1237 → 1250.

> Take-home: the path-network's `shuffle(rng, targets)` step
> (FW2011...md:156) is FW-prescribed but breaks the analytic
> `u(z̄) = ū(z)` invariant for real-coef real-IC ODEs.  The opt-in
> kwarg restores bit-exact Schwarz symmetry by walking only upper-
> half + on-axis canonical representatives and mirroring lower-half
> via `conj`.  Roughly halves wall-time on full-plane fills.

## What changed

  - `src/PathNetwork.jl`:
    - New helper `_solve_with_schwarz_reflection` (~70 LOC) — validates
      IC-on-real-axis (`Im(z₀) = Im(u₀) = Im(u'₀) = 0` to within
      `10*eps(T)`), canonicalizes input cells via `complex(real(z),
      abs(imag(z)))` (collapses ±0.0 ambiguity), recurses with
      `enforce_real_axis_symmetry=false` on the unique upper targets,
      builds output `grid_u` / `grid_up` by indexing-and-conjugating.
    - `path_network_solve` gains the kwarg.  Early-return dispatches
      to the helper when `true`; default `false` preserves the FW 2011
      algorithm verbatim.
    - Module-header docstring gains a new chapter "Schwarz-reflection
      symmetry (opt-in)" between "Per-node canonical Padé" and
      "Tier-3-to-5 deferrals".

  - `test/pathnetwork_test.jl` — 4 new testsets:
    - **PN.6.1** — bit-exact `u(z̄) = ū(z)` on a 5×5 conjugate-symmetric
      grid; checks `max_asym_u == 0.0` AND `max_asym_up == 0.0` over
      20 off-axis pairs (on-axis pairs are trivial and skipped).
    - **PN.6.2** — visited tree confined to `Im(z) ≥ 0` on a grid
      forcing multi-step walks (off-axis targets beyond `h`).  Also
      checks `length(visited_z) > 1` to ensure walking actually
      happened (catches mutations that disable the filter without
      breaking trivial-coverage grids).
    - **PN.6.3** — input grid order preserved in `sol.grid_z`.
    - **PN.6.4** — fail-fast on off-axis IC (`zspan[1]`, `u₀`, `u'₀`).

  - `examples/tritronquee_3d.jl` — drop the manual upper-half walk +
    conjugate mirror; pass the full grid through with
    `enforce_real_axis_symmetry=true`.  Net `-15 LOC` user-space code,
    same PNG output.

## Mutation-proof procedure

Two mutations applied + reverted; both bit as predicted.

**Mutation G** — in the output-build loop, drop `conj()` in the
`imag(z) < 0` branch so lower-half cells inherit the upper-walk value
as-is.  Verified bite: PN.6.1 lines 257-258 RED (`max_asym_u` /
`max_asym_up` non-zero).  PN.6.2-4 stayed GREEN (they don't depend on
the mirror itself).  Confirms the conj mirror is load-bearing for the
Schwarz invariant.

**Mutation I** — in `_solve_with_schwarz_reflection`, replace
`complex(real(z), abs(imag(z)))` with `complex(real(z), imag(z))` (drop
`abs`).  The upper-canon collapse breaks: lower-half cells map to
themselves, the recursive walk visits both halves, and the conj branch
applies to an asymmetric walk's value.  Verified bite: PN.6.1
lines 257-258 RED AND PN.6.2 line 279 RED.  Validates BOTH the symmetry
mirror AND the upper-half filter as load-bearing.

Both restored before commit; GREEN at 1250 / 1250.

## Design decisions

**Canonical-representative collapse via `abs(imag(z))`.**  An input grid
may contain `0.5 + 0.0im` and `0.5 - 0.0im` as distinct entries (Julia's
`isequal(+0.0, -0.0)` is `false`); the dict-lookup in the mirror loop
would then KeyError on the second one.  `complex(real(z), abs(imag(z)))`
collapses both to `0.5 + 0.0im` (positive zero) before unique/dict-
keying, eliminating the signed-zero footgun.

**Validation policy.**  The kwarg's correctness depends on TWO
preconditions: (a) real coefficients in `f`, and (b) real IC on the
real axis.  (a) cannot be checked at the path-network layer without
symbolic analysis of `f` — the kwarg pushes that burden onto the caller
(documented in the docstring + module chapter).  (b) IS checked
(`abs(imag(prob.zspan[1])) ≤ 10*eps`, same for `y0[1]` and `y0[2]`) and
throws `ArgumentError` per CLAUDE.md Rule 1.  PN.6.4 covers all three
checks.

**Recursive call into `path_network_solve(...; enforce_real_axis_symmetry=false)`.**
Rather than duplicating the Stage-1/Stage-2 loop with a filtered grid,
the helper recursively calls the existing implementation on the upper-
canon targets.  Clean, no code duplication, and any future changes to
the core algorithm automatically inherit.  The recursion terminates
in one step (the recursive call has `false` explicitly).

**On-axis output handling.**  On-axis input cells receive the upper-
walk's complex value as-is (no `real()` projection).  The path-network
walks complex paths to reach real-axis targets and may produce small
non-zero imag noise; we don't enforce `Im(u) = 0` exactly on the real
axis.  Callers who need strictly-real on-axis values can take `real(...)`
themselves.  Document this as the kwarg's contract: the bit-exact
Schwarz invariant holds for `Im(z) ≠ 0` pairs; on-axis cells inherit
the walker's natural output.

## Frictions surfaced

  - **Signed-zero KeyError in the test loop** (caught on first GREEN
    attempt).  Tests built `idx_of = Dict(z => i)` keyed by `sol.grid_z`
    (positive zeros from `range(...) + im*y`) and looked up `conj(z)`
    (which gives `-0.0im` on the on-axis row).  `isequal(+0.0im, -0.0im)
    == false` ⇒ KeyError.  Fix: skip on-axis cells in the pair check
    (trivially `z == conj(z)`).  Also recorded as the impl design point
    above.

  - **PN.6.2's first grid was too cheap to bite Mutation I.**  Initial
    grid had all targets within `h=0.5` of IC at `z=0` ⇒ 0 steps walked
    ⇒ `visited_z == [0+0im]` regardless of mutation.  Rewrote to use
    off-axis targets > `h` from IC (`0.7-1.2im`, etc.) so Stage-1
    actually grows the visited tree.  Lesson: when designing tests
    that assert properties of the visited tree, the grid must force
    walking.

## Beads

  - `padetaylor-dtj` — closed in this session.

## Pointers

  - `src/PathNetwork.jl` lines 60-89 (module-header chapter), lines
    132-150 (function docstring + signature), lines 274-345 (helper).
  - `test/pathnetwork_test.jl` PN.6.1-6.4 + mutation-proof block at
    end-of-file.
  - `examples/tritronquee_3d.jl` — the user-space workaround retired
    in favor of the kwarg.
  - `docs/worklog/014-pathnetwork-symmetry-debug.md` — the original
    bug discovery + user-space workaround that motivated this work.
  - `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md`
    line 156 (FW's shuffle prescription, preserved as default behavior).
