# Worklog 081 — Diagnostics silent sample loss

**Date:** 2026-09-29 (pace session, worktree `w1-diag`, branch `pace/w1-diag`).
**Scope:** `padetaylor-orb5` and the mask/accounting half of `padetaylor-qdsm`.
**Base:** `b9f4132`; changes intentionally uncommitted for the orchestrator.

## Verified ground truth and stale pointers

Read `CLAUDE.md` and `AGENTS.md` before editing. At the base commit the
beads' principal pointers are accurate: the old mask is
`ext/PadeTaylorDiagnosticsExt.jl:153-161`, consumed at :226-227, and its
selected node cloud is triangulated at :239-241. The fallback at :157
is the implicit ζ strip. The exceptions are silently skipped at :276-285,
with nonfinite values skipped at :286-287 (the bead's :277-286 slightly
understates the latter span). The four aggregates are :312-315.
`src/Diagnostics.jl:152-168` is the old report, :174-181 its public mask
description, and :50-58 the module-docstring bullet, not implementation.
`test/diagnose_test.jl:111-124` is DG.5 and remains untouched.
`docs/worklog/077-fig47-path-dependence-seam-rootcause.md:188-190` is the
original next-step paragraph; it now names qdsm and this implementation.
These old positions are explicitly base-commit citations, not current ones.

FFW's strip applies to the map `z=exp(ζ/2)`:
`references/markdown/FFW2017_painleve_riemann_surfaces_preprint/FFW2017_painleve_riemann_surfaces_preprint.md:101-103`.
The original probe prints its selected population
(`external/probes/loop-closure-fig1/probe.jl:152-163`); diagnostics omitted
that protection. The scalar evaluator throws a pole `DomainError` at
`src/PadeStepper.jl:376-380`. Its low-to-high coefficient convention is
`src/RobustPade.jl:120-130`; read `external/chebfun/padeapprox.m:1-28` too.
FW's 20×20 independent z-plane windows are documented at
`references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:147`.

## Changes and why

Choose qdsm option (b): branchless metadata retains every node, rather
than guessing a frame (`ext/diagnostics_scalar.jl:10-20`). Explicit
branched sheet-0 filtering and DG.5's nonzero-sheet rejection remain.
`min_retained_fraction=1.0` rejects any unrequested node loss; a caller
can explicitly budget a known subset (:22-35, :44-66). No new mask or
sheet keyword competes with `padetaylor-s1q`.

Report counts and display expose nodes retained/dropped, non-tree edges
candidate/evaluated/dropped, pole exceptions and nonfinite results
(`src/diagnostics_reports.jl:57-89`, :94-108). Every failed edge counts
once, including failures at either endpoint. Unknown exceptions propagate;
all-candidate failure throws with counts and a suggestion
(`ext/diagnostics_edges.jl:25-40`, :51-69).

Off-sheet losses in the current algorithm are parent links excluded
before tree subtraction, not midpoint failures. Count these as
`n_tree_edges_off_sheet` (`ext/diagnostics_geometry.jl:80-92`) separately
from the non-tree candidate denominator. No fictitious all-sheet Delaunay
population is claimed. This preserves the original sheet-local graph.

Scalar/vector methods share geometry, accounting and aggregation; vector
coverage remains complete (`ext/diagnostics_vector_adapter.jl:45-63`).
The internal `allow_degenerate` flag preserves the vector adapter's
existing <3-node empty report and the scalar method's existing
Delaunay insufficient-point error; the refactor does not lift that
scalar limitation. If scalar one-node diagnostics is required, the
orchestrator should file a separate bead for its documented error policy.
Depth/LCA/category helper bodies were compared against `git show HEAD`:
exactly preserved. All touched Julia code and the new test are below 200
physical lines, including chapter documentation. ADR-0016 links its
sampling amendment; README documents the public coverage budget and counts.

## Measurements and independent oracle

| Fixture | Measured result |
| --- | --- |
| Im z ∈ [-50,50], x ∈ {-1,0,1}: 303 nodes | Old strip 39/303 = 12.871287%; 264/303 = 87.128713% dropped; corrected mask 303/303 |
| Explicit sheet metadata with 3/4 principal nodes | Default rejects; budget 0.75 accepts 3 retained / 1 dropped; budget 0.8 rejects; 1 off-sheet parent link |
| Non-cocircular five-node constant-field fixture | 7 Delaunay edges, 2 Delaunay tree edges, exactly 5 candidate non-tree edges |
| One exact midpoint denominator root, at either endpoint | 1 thrown / 0 nonfinite / 1 dropped / 4 evaluated |
| NaN patch at node 2 / node 5 | 3 / 2 nonfinite edges; 2 / 3 evaluated survivors with max ΔP_rel = 0 exactly |
| Two-component vector constant field, NaN at node 2 | 5 nodes retained, 5 candidates, 3 nonfinite dropped, 2 evaluated; max ΔP_rel = 0 exactly |
| Triangle with its only non-tree midpoint on a pole | 1 candidate, 1 thrown; ArgumentError with Suggestion |

The field oracle is exactly `u(z)=1`, or vector `(1,2)`; its finite
midpoint disagreement is zero, independently of diagnostics code.
For points `(0,0),(4,0),(3,4),(0,2),(5,1)`, independently enumerate all
triples in Python using `fractions.Fraction` circumcentres, retain those
with strictly empty circumcircles, and confirm no fourth point lies on a
circle. The three triangles are `(1,2,4),(2,3,4),(2,3,5)`. Their union
gives the seven edges; removing root-star links leaves
`{(2,3),(2,4),(2,5),(3,4),(3,5)}`. The test pins this known set at
`test/diagnostics_sampling_test.jl:80-84`; graph degree then gives exact
nonfinite losses. Initial cocircular fixture failed 5 assertions because
valid triangulations differed; replacing its geometry fixes that cause,
without changing a tolerance or deleting assertions.

## Test execution and mutation evidence

Julia 1.12.3, one process at a time, private depot. Root standalone
loading initially failed because Delaunay is an extra/weak dependency.
It loads offline from the existing Fig 1 probe environment by stacking:

```bash
export JULIA_LOAD_PATH="$PWD:$PWD/external/probes/loop-closure-fig1:@stdlib"
julia() { /home/tobiasosborne/.julia/juliaup/julia-1.12.3+0.x64.linux.gnu/bin/julia "$@"; }
JULIA_DEPOT_PATH="$PWD/.depot:$HOME/.julia" julia --project=. test/diagnostics_sampling_test.jl
```

The shell function bypasses juliaup's attempted read-only configuration
lock; it uses the already installed Julia binary. No dependency install,
resolution, network, full suite or shared-depot writes were performed.
Mutation logs are local ignored `.depot/diagnostics-mutation-*.log` files.
After final restoration the new suite passed **62 / 0 fail / 0 error / 0 broken**
in 17.1 s. Each mutation was applied alone, tested with the same command,
and restored before the next one. All three mutated files byte-match their
saved correct snapshots (`cmp` exited 0 for each). The complete procedure
is also recorded in `test/diagnostics_sampling_test.jl:153`.

| Mutation | Exact perturbation | RED output; pass/fail/error |
| --- | --- | --- |
| A: mask | `ext/diagnostics_scalar.jl:13`: `trues(length(visited_z))` → `[(-2π < imag(z) ≤ 2π) for z in visited_z]` | `Evaluated: 39 == 303` twice; `retained 39 / 303 nodes (12.871287128712872%)`; 52/2/1 |
| B: node-loss guard | :28: prepend `false &&` to the `if` condition | `Expression: quality_diagnose(sol)` / `No exception thrown`; 50/1/1 |
| C: pole count | `ext/diagnostics_edges.jl:29`: `n_threw += 1` → `n_threw += 0` | `Expression: pole.n_edges_threw == 1` / `Evaluated: 0 == 1`; 55/5/1 |
| D: nonfinite count | :34 and :39: `n_nonfinite += 1` → `n_nonfinite += 0` | `bad.n_edges_nonfinite == length(incident)` / `0 == 3` and `0 == 2`; vector `0 == 3`; 52/10/0 |
| E: off-sheet links | `ext/diagnostics_geometry.jl:87`: `n_off_sheet += 1` → `n_off_sheet += 0` | `r.n_tree_edges_off_sheet == 1` / `0 == 1`; 61/1/0 |
| F: broad catch | `ext/diagnostics_edges.jl:28`: `err isa DomainError || rethrow()` → `true` | `Expression: run()` / `No exception thrown`; 60/1/1 |
| G: total failure | :55: prepend `false &&` to the `if` condition | `Expression: quality_diagnose(sol)` / `No exception thrown`; 59/1/1 |
| H: invalid fraction | `ext/diagnostics_scalar.jl:23`: validation expression → `true ||` | invalid-fraction calls: `No exception thrown` twice; 60/2/0 |

Every mutant exited 1 and had 0 broken tests. Error entries are the
expected-error capture helper loudly aborting after a missing exception.
The off-sheet mutation was repeated after preserving the scalar/vector
degenerate-input distinction: again 61/1/0 RED (19.5 s), restored, then
62/0/0/0 GREEN. The final restored source retains all eight mechanisms.
No localized-gate mutation is claimed, because that gate is not implemented.

Final standalone checks (all serial; no full/fast aggregate gate invoked):

| File | Pass | Fail | Error | Broken | Time |
| --- | ---: | ---: | ---: | ---: | --- |
| `test/diagnostics_sampling_test.jl` | 62 | 0 | 0 | 0 | 17.1 s |
| `test/diagnose_test.jl` | 52 (32 scalar + 20 vector) | 0 | 0 | 0 | 20.3 s + 5.5 s |
| `test/pathnetwork_test.jl` | 123 | 0 | 0 | 0 | 1m18.0s |

`git diff --check` passes. Read-only verification confirms every forbidden
file is unchanged and every touched Julia file is ≤175 physical lines.
Register the new file above in `test/runtests.jl` during orchestration;
that file was deliberately not edited here. Global expected-broken count
was not measured by this permitted subset of test files.

## Deferred qdsm gate design and orchestration follow-ups

The second half is not implemented or validated on a Fig 4.7 solve.
Keep `padetaylor-qdsm` open until a full-node seam-versus-p99 fixture
and a mutation replacing the gate with p99 are measured. The trigger is
using diagnostics as a pass/fail criterion for Fig 4.7 or refinement.
Worklog 077 :151-175 motivates field/path-independence checks; :81-84
records 3523/10201 cells differing by >1% between seeds.

Proposed design: bin non-tree midpoint coordinates on an explicitly
specified spatial grid, report each bin's maximum ΔP_rel, evaluated and
dropped counts, and worst edge/locus. Separate in-disc (`extrap_max ≤ 1`)
agreement from extrapolated comparisons. Gate on the worst eligible bin
at a caller-specified tolerance; coverage failures cannot pass the gate.
Keep the four existing aggregates as context. Include a known localized
seam affecting <1% of comparisons, where p99 stays below threshold but
the appropriate local maximum exceeds it, then verify on the actual
full-node Fig 4.7 field and two seeds. Do not claim that spatial binning
alone adds sensitivity: maximum over bin maxima equals the existing
global maximum if all edges are eligible. The new work must establish
a useful pass/fail contract, spatial locus and trustworthy population.
Moreover, locally agreeing patches can still be wrong together; the
two-seed field invariant is necessary independent evidence, not something
a single-walk midpoint maximum proves.

The orchestrator must register `test/diagnostics_sampling_test.jl`, link
`padetaylor-6pj`, `padetaylor-adn`, and `padetaylor-s1q` to qdsm, record orb5
completion and qdsm's remaining gate scope, and run permitted aggregate
gates before committing. Tracker files, `bd`, commits, pushes, full-suite
commands, `HANDOFF.md`, `CHANGELOG.md`, and `test/runtests.jl` were excluded
by the common brief. The positional DiagnosticReport constructor gained
eight fields; downstream manual constructors require adjustment.
