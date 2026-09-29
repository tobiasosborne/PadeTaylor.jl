# ADR-0016 Amendment — Explicit diagnostics sample coverage

**Status:** Accepted (2026-09-29), mask/accounting portion only.
**Beads:** `padetaylor-qdsm` (first half), `padetaylor-orb5`.
**Parent:** [ADR-0016](0016-diagnostics-extension.md).

## Ground truth

FFW 2017 defines the strip `-2π+4πs < Im ζ ≤ 2π+4πs` for the specific
map `z=exp(ζ/2)`, not for arbitrary z-plane nodes:
`references/markdown/FFW2017_painleve_riemann_surfaces_preprint/FFW2017_painleve_riemann_surfaces_preprint.md:101-103`.
The original Fig 1 probe used that frame and explicitly printed its
retained population (`external/probes/loop-closure-fig1/probe.jl:152-163`).
The diagnostics extension copied the predicate but omitted the accounting.
FW 2011 uses independent 20×20 z-plane windows for Fig 4.7
(`references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:147`).
Empty sheet metadata therefore cannot justify a coordinate strip.

## Decision

Choose qdsm option (b): branchless walks retain every node. Do not infer
the frame from empty `visited_sheet`, and do not introduce `sheet_mask`.
Branched walks keep the existing `[0]` predicate; the `sheet=0` API and
`sheet != 0` error remain. Multi-sheet/cut-aware triangulation remains
`padetaylor-8py` work; the sheet-tuple API is still `padetaylor-s1q` work.
The old probe baseline must be compared on an explicitly selected
ζ-frame node cloud, not by reinstating an implicit strip in diagnostics.

The scalar method adds `min_retained_fraction::Real=1.0`. Any node loss
is rejected by default. A caller intentionally diagnosing a sheet-0
subset can supply its minimum fraction explicitly. The fraction must be
finite and in [0,1], and errors include observed counts plus a suggestion.
This budget controls coverage only; it does not permit other sheets.

`DiagnosticReport` gains these exact population fields:

| Field | Population |
| --- | --- |
| `n_nodes` | All input visited nodes |
| `n_nodes_retained` / `n_nodes_dropped` | Nodes selected / excluded by the mask |
| `n_edges_candidate` | Non-tree Delaunay edges on the retained cloud |
| `n_edges` | Successfully evaluated candidate edges (unchanged meaning) |
| `n_edges_dropped` | Candidate edges that could not be evaluated |
| `n_edges_threw` | Candidate edges raising a Padé pole `DomainError` |
| `n_edges_nonfinite` | Nonfinite endpoint values or disagreement statistics |
| `n_tree_edges_off_sheet` | Parent links excluded because either endpoint was masked |

Each failed candidate counts once, including failure at either endpoint.
Unknown exceptions propagate instead of being silently classified as poles.
All-candidate evaluation failure throws with counts and a suggestion.
A valid empty candidate population keeps the pre-existing empty report
convention; fewer than three vector nodes keep that convention, while
scalar insufficient-point errors remain unchanged. Quantiles and category
counts remain over successful edges.
The text display prints sample coverage before quality statistics.

Off-sheet parent links are counted separately: they are tree edges, not
failed members of the non-tree evaluation population. A triangulation
only on retained nodes cannot count non-tree edges on excluded sheets
without inventing another graph. We keep the original sheet-local
triangulation; no cross-sheet edge count or cut detection is claimed.

Both scalar and vector methods share geometry, evaluation accounting and
aggregation. The vector method retains all nodes and takes no fraction
or sheet keyword. Helpers are split into files below 200 physical lines.
Adding report fields changes direct positional construction of the public
struct; repository callers use the updated extension constructors.

## Validation and deferred gate

`test/diagnostics_sampling_test.jl` uses the independent constant-field
oracle `u(z)=1`, known midpoint roots, and a 303-node fixture spanning
Im z ∈ [-50,50]. The old strip retains 39/303 (12.8713%); the corrected
mask retains 303/303. Coverage rejection, pole/nonfinite accounting at
both endpoints, the vector adapter and all-candidate failure are tested.
Mutation evidence and exact test results are in
`docs/worklog/081-diagnostics-silent-drops.md`.

The localized maximum gate is **not shipped** in this amendment.
`padetaylor-qdsm` remains open for that second half; retain median/p90/p99/max
as context. See worklog 081 for a spatial-bin design and its required
Fig 4.7 seam-versus-p99 regression. `padetaylor-6pj` and `padetaylor-adn`
can now measure full z-frame coverage but must not treat these aggregates
as a validated seam acceptance gate. Tracker changes belong to the pace
orchestrator, since this engineer is forbidden to edit `.beads/`.
