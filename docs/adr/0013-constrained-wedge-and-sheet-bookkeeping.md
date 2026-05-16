# ADR-0013 — Constrained-wedge routing + per-branch sheet bookkeeping (FFW 2017 §2.2.2)

**Status**: Accepted (2026-05-16) | **Bead**: `padetaylor-bho` | **Worklog**: 042

## Context

The path-network walker shipped in ADR-0004 (Tier 2) advances by
selecting the wedge candidate that minimises `|u|` (or steepest-descent).
There is no notion of a **branch cut**: the walker is free to choose any
direction whose Padé step happens to be valid.

For PVI in the ζ-plane (FFW 2017 §2.2 / SheetTracker primitives in
worklog 018) — and for any multi-valued solution computed in a frame
where the equation's fixed branch points sit at finite locations —
the walker MUST respect the Riemann-sheet structure of the analytic
solution:

  - On a chosen sheet, the path may not cross a branch cut (the
    cut is the artefact of the principal-branch convention; crossing
    it would land on a different sheet without the bookkeeping ever
    noticing).
  - Crossing a cut DELIBERATELY is the way the walker moves between
    sheets — and is exactly how FFW Fig 2 column 3 / Fig 3 are
    rendered. The walker counts these crossings to label each visited
    node with its sheet index.

The FFW 2017 §2.2.2 algorithm (`references/markdown/
FFW2017_painleve_riemann_surfaces_preprint/
FFW2017_painleve_riemann_surfaces_preprint.md:163-189`) describes
both modes:

> "To move onto the sheet in the top-right frame the paths move
> through the branch cut by running only in a counterclockwise
> direction around `ζ = 0`, as indicated by the arrows. The paths
> then run in both directions around the branch point at `ζ = 2πi`.
> ... none of the paths overstep the branch cuts, indicated by solid
> lines, on a given sheet." — FFW md:178

The v1 SheetTracker primitives (`winding_delta`,
`accumulate_winding`, `sheet_index`; worklog 018) compute sheet
indices AFTER a path-network walk has completed — they cannot prevent
overstepping during the walk, nor can they thread sheet bookkeeping
into the wedge selector. Bead `padetaylor-bho` materialises the
production-grade version: a `path_network_solve` extension whose
wedge selector consults a cut-crossing predicate and (in the opt-in
`:cross_branch` mode) maintains per-branch sheet counters as nodes
are appended to the visited tree.

This is **step A4** of the 11-step FFW 2017 figure-reproduction
plan (steps A1, A2, A3, A6 shipped; B1, B2, B3, B4 shipped). A4 plus
A5 (sheet-aware Stage-2) together unblock B5 (FFW Figs 2, 3, 7 — the
multi-sheet PVI rendering).

## Decision

Ship **two opt-in kwargs** on `PathNetwork.path_network_solve`:

```julia
path_network_solve(prob, grid;
    branch_points     :: Tuple{Vararg{Complex}}        = (),
    branch_cut_angles :: Union{Real, Tuple{Vararg{Real}}} = π,
    cross_branch      :: Bool                            = false,
    # ... existing kwargs unchanged ...
)
```

plus a **new helper module** `src/BranchTracker.jl` that owns the
cut-crossing predicate, per-step winding-delta computation, and the
sheet-counter accumulator. The helper module ships ~150 LOC and
keeps `path_network_solve`'s Stage-1 inner loop's growth to ~50 LOC
of branch-aware threading.

**Algorithm spec** (FFW md:178-188 verbatim):

  - **Branch cut model**: each branch point `b_k ∈ ℂ` has a
    half-line cut `{b_k + t · exp(i · α_k) : t ≥ 0}`, where `α_k`
    is `branch_cut_angles[k]` (or the scalar `branch_cut_angles`
    broadcast to all branches). Default `α_k = π` matches Julia's
    `log` convention: cuts run from each branch point along the
    negative-real direction.

  - **Cut-crossing predicate** (`BranchTracker.segment_crosses_cut`):
    given segment `S = [z_cur, z_new]` and cut `C = [b, b + ∞·e^{iα}]`,
    solve the standard 2×2 parametric intersection
    `z_cur + t(z_new - z_cur) = b + s·e^{iα}` for `(t, s) ∈ ℝ²`,
    and return `true` iff `t ∈ (0, 1)` and `s > 0` (open interval
    on `t` so steps that *land on* the cut endpoint are not
    counted; closed on `s` similarly to avoid double-counting
    the branch point itself). Parallel segments
    (determinant zero) return `false` — degenerate.

  - **Refuse mode (`cross_branch = false`, default)**: in the
    wedge loop, after each candidate `z_new` is computed, check
    every branch cut via `segment_crosses_cut`. If ANY cut is
    crossed, mark the candidate as forbidden (treat like a
    failure). If ALL 5 wedge candidates are forbidden, throw an
    `ErrorException` whose message points the caller at
    `cross_branch = true` (per CLAUDE.md Rule 1). Visited-sheet
    counters stay at the root's value (`(0, 0, ...)`) for every
    node — no traversal accumulates winding.

  - **Cross mode (`cross_branch = true`)**: candidates are
    allowed to cross cuts. Per step, for each branch `k`, compute
    the winding delta
    `Δθ_k = wrap((arg(z_new - b_k) - arg(z_cur - b_k)) + π) - π`
    (normalised to `(-π, π]`; the SheetTracker `winding_delta`
    primitive shipped in worklog 018 already does this). If the
    step crosses cut `k` (per the same predicate), the sheet
    counter for `b_k` is updated by `sign(Δθ_k)`:
    counter-clockwise crossing (`Δθ_k > 0`) bumps `+1`,
    clockwise (`Δθ_k < 0`) bumps `-1`. Stored at
    `visited_sheet[node_idx][k]`. Per-step caller contract: step
    magnitude must be small enough that `|Δθ_k| < π` at every
    branch (otherwise winding becomes ambiguous; same precondition
    as the SheetTracker primitives).

  - **Storage**: `PathNetworkSolution` grows a new field
    `visited_sheet :: Vector{Vector{Int}}` of length
    `length(visited_z)`, each inner vector of length
    `length(branch_points)`. Empty inner vectors (`Int[]`) when
    `branch_points = ()`, preserving wire-format compatibility
    with existing tests. The root node's entry is the
    `initial_sheet` kwarg (default zeros).

**API summary** (the public surface):

```julia
sol = path_network_solve(prob, grid;
    branch_points     = (0im, 2π*im, -2π*im),
    branch_cut_angles = π,                          # all 3 cuts along arg = π
    cross_branch      = true,
    # ... other kwargs ...
)

sol.visited_sheet[k] :: Vector{Int}    # length 3; per-branch sheet
                                       # index at node k
```

Backward compat: with `branch_points = ()` (the default), the new
code paths are byte-equivalent to the pre-A4 walker — no cut checks,
no sheet bookkeeping, no `visited_sheet` allocation cost (an empty
`Vector{Vector{Int}}` is still allocated for shape consistency, ~16
bytes). All 2262 existing tests must remain GREEN. This is the
load-bearing invariant.

## Why these specific choices

### Why half-line cuts not arbitrary segments?

Three reasons:

1. **Matches FFW practice**: every FFW 2017 figure that exhibits
   cut-respecting routing (Figs 2, 3, 7) operates in the ζ-plane,
   where the natural cuts are half-lines from each branch lattice
   point (FFW md:178, md:185-189 sheet parametrisation).

2. **Matches Julia convention**: `log(z)` cuts along `arg = π`;
   `√z` cuts along `arg = π`. Users reaching for branch cuts in
   numerical code already think in this idiom.

3. **Per-branch direction lets the model cover finite cuts**:
   when finite cuts (e.g., z-plane PVI's `z ∈ (0, 1)`) become
   necessary, the user can supply two anti-parallel half-line cuts
   (`α_0 = 0` from `b_0 = 0`, `α_1 = π` from `b_1 = 1`) — the
   geometric intersection coincides with the segment. This is a
   v2 ergonomic-only lift if it matters.

Arbitrary line-segment cuts (the second option considered) would
add API surface for no v1 benefit; user-predicate `forbidden_step`
(the third) would throw away the structural branch-point info
needed for sheet bookkeeping.

### Why a new helper module rather than inlining into PathNetwork.jl?

PathNetwork.jl is already at 744 LOC against the Rule 6 cap of 200
(open bead `padetaylor-qum`). Adding ~150 LOC of branch-tracking
inline pushes it to ~900. A new `src/BranchTracker.jl` module:

  - Keeps PathNetwork.jl's growth to ~50 LOC (kwarg threading +
    one helper call per wedge candidate).
  - Owns the geometric primitives in one cohesive place (cut
    intersection, winding-delta-with-cut-detection, sheet-counter
    update), independently testable.
  - Does NOT delete the existing SheetTracker primitives (which
    operate on already-walked paths); the two modules occupy
    complementary niches — BranchTracker is the **walker-side**
    enforcement layer, SheetTracker is the **caller-side**
    analysis layer.

This is the same split discipline as ADR-0004 (PathNetwork sibling
to Problems) and ADR-0014 (IVPBVPHybrid composes PathNetwork +
BVP). Composition over inheritance.

### Why per-branch sheet TUPLE rather than scalar?

FFW md:188 parametrises PVI sheets by `(θ_k, φ_k)` — TWO
independent integer indices, one for each of the two branch points
(`z = 0` and `z = 1`). The tuple representation captures this
structurally: a walker that circumambulates only one branch
increments only that index. A scalar "sheet number" would have to
fuse the two via some convention (e.g., `2·k + l` for PVI's two
branches) and would not generalise to PIII's lattice of
`ζ = 2π·i·k` branches.

The cost is one extra integer per visited node per branch (8 bytes
each at Int64). For a typical FFW Fig 2 walk with ~3000 visited
nodes and 1 branch in scope, that's 24 KB total — negligible.

### Why fail-loud when all 5 wedges forbidden?

Per CLAUDE.md Rule 1: crashes-with-context beat truthful-looking
lies. The natural alternatives are silently returning the best
candidate even if forbidden (would corrupt sheet bookkeeping), or
expanding the wedge to find an allowed direction (would mask user
error). The error message names `cross_branch = true` as the most
common fix.

## Consequences

### Test coverage

A new test file `test/branch_tracker_test.jl` covers the helper
module primitives (segment-cut intersection on canonical cases,
winding-delta accumulator, sheet-counter update direction). A
second new file `test/path_network_branch_test.jl` covers the
walker-level integration (refuse mode rejection, cross mode
counter increment, backward-compat byte-equality on the default
`branch_points = ()` path, fail-loud message content).

Target: ~80 new GREEN assertions; three mutations (a flip on the
cross-direction sign, a drop of the cut check in refuse mode, a
zero-out of the per-branch winding accumulator) must each bite a
distinct subset of the new tests.

### Documentation

  - This ADR-0013.
  - Worklog 042 documenting the shipping session.
  - `path_network_solve`'s docstring grows three kwarg paragraphs
    (already long; this pushes Rule 10 / literate-programming
    discipline toward the LOC-cap split deferred to bead
    `padetaylor-qum`).
  - `figure_catalogue.md §6` T5 row updated from PARTIAL to a
    further intermediate status pending A5 + B5.
  - HANDOFF.md refresh at session close.

### Open follow-ups (v2 work, not in this bead)

  - **Finite cut segments**: ergonomic kwarg
    `branch_cut_segments::Vector{Tuple{Complex,Complex}}` for
    callers who want to specify a finite cut endpoint pair directly
    rather than expressing it as two anti-parallel half-lines. File
    a P3 bead if a downstream caller asks for it.

  - **A5 sheet-aware Stage-2**: separate bead. The Stage-2
    barycentric lookup currently picks the nearest visited node by
    Euclidean distance; A5 restricts the pool to visited nodes
    whose `visited_sheet[k]` matches the query's sheet. This is
    necessary for FFW Fig 2/3 multi-sheet rendering — A4 alone
    populates `visited_sheet` but A5 is what makes the lookup
    respect it.

  - **B5 figures (FFW Figs 2/3/7)**: separate beads, gated on
    A5. Each figure script will exercise the A4+A5 stack on its
    own parameter set.

  - **PathNetwork.jl LOC-cap split** (`padetaylor-qum`): the
    growing file should split into `Stage1.jl` / `Stage2.jl` /
    `Symmetry.jl` / coordinator. A4 makes the cap pressure worse
    but the split itself is a separate refactor bead.

  - **`path_network_solve(::PainleveProblem)` forwarding** for
    `branch_points` / `cross_branch`: currently the Painleve
    layer's forwarding methods don't expose these kwargs. Trivial
    to add when a B5 figure-script calls them.

## References

  - **FFW 2017 §2.2.2** — `references/markdown/FFW2017_painleve_riemann_surfaces_preprint/FFW2017_painleve_riemann_surfaces_preprint.md:163-189` (circumambulation algorithm), md:178 (the "none of the paths overstep" specification + the deliberate-crossing description), md:185-189 (sheet parametrisation via `(θ_k, φ_k)` index ranges).
  - **Worklog 018** — `docs/worklog/018-sheet-tracker-pVI.md` §"What is NOT shipped" — the v1 deferral that this ADR closes.
  - **Worklog 041 / ADR-less A3** — η-plane PVI RHS; the immediately preceding A4 substrate.
  - **ADR-0004** — `docs/adr/0004-path-network-architecture.md` (PathNetwork as sibling to Problems; the composition discipline that motivates the helper-module split).
  - **ADR-0011 / ADR-0012** — `docs/adr/0011-adaptive-pade-step.md`, `docs/adr/0012-non-uniform-stage-1-nodes.md` — the two preceding opt-in PathNetwork kwargs whose orthogonality pattern A4 follows.
  - **`src/SheetTracker.jl`** — `winding_delta`, `accumulate_winding`, `sheet_index` (worklog 018) — the after-the-fact primitives that BranchTracker complements with walker-side enforcement.
  - **`src/PathNetwork.jl`** — the file the new kwargs land in; bead `padetaylor-qum` tracks the orthogonal LOC-cap split this work makes more pressing.
