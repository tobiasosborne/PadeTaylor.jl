# ADR-0012 — Non-uniform Stage-1 node placement (FFW 2017 §2.1.2)

**Status**: Accepted (2026-05-15) | **Bead**: `padetaylor-1a3` | **Worklog**: 035

## Context

The FW 2011 path-network walker advances by a **constant** step magnitude
`h`. ADR-0011 retired the `h ≡ const` floor on the *temporal* axis with
the `:adaptive_ffw` controller, which shrinks `h` adaptively under a
per-step truncation-error estimate. The *spatial* axis — where the
walker chooses to land — remained governed by the same constant
magnitude.

FFW 2017 §2.1.2 (`references/markdown/FFW2017_painleve_riemann_surfaces_preprint/
FFW2017_painleve_riemann_surfaces_preprint.md:67-72`) observes that
this is catastrophic for `P̃_III` / `P̃_V` solutions whose pole densities
**grow exponentially** with `Re ζ` under the standard exponential
transformations (CoordTransforms shipped in worklog 017):

> "the exponential transformations used to arrive at `P̃_III` and
> `P̃_V`… [imply that] the pole density will increase rapidly on the
> region Re ζ ≫ 0. For simplicity we have therefore chosen Stage 1 node
> sets with a node separation function R(ζ) that decreases linearly with
> Re ζ" — FFW md:72.

For their headline PIII solution (md:101 Figure 1, `(α,β,γ,δ) =
(-1/2,-1/2,1,-1)`, `(z,u,u') = (1, 1/4, 1)`), FFW prescribe
`R(ζ) = (8 - Re ζ)/20`. A uniform-`h` walker spends most of its budget
on the smooth low-density region and arrives under-resolved at the
high-density region.

Worklog 017 §"What is NOT shipped" deferred this with two options:
**(a)** a higher-level grid-generator that wraps `path_network_solve`,
or **(b)** a kwarg on `path_network_solve` itself. Bead `padetaylor-1a3`
materialises this deferral.

This is **step A2** of an 11-step plan to reproduce FFW 2017's seven
figures. Step A1 (adaptive `h`, ADR-0011) shipped in commit `7c7cac7`;
this ADR stacks orthogonally on that work.

## Decision

Ship the kwarg-on-walker (option **b**): an opt-in `node_separation ::
Union{Nothing, Function} = nothing` kwarg on
`PathNetwork.path_network_solve`. When `nothing` (the default),
behaviour is byte-identical to the FW 2011 / ADR-0004 / ADR-0011
combined stack — the load-bearing backward-compat invariant. When a
function `R(z) -> real`, the per-step walker step magnitude is `R(z_cur)`
at every Stage-1 inner step.

**Algorithm spec** (verbatim FFW md:67-97):

  - At each step of a per-target walk, before the wedge candidate
    evaluation, the initial step magnitude is `h_seed = R(z_cur)`
    rather than the controller's memory seed (`h_cur`).
  - Under `step_size_policy = :fixed`, `h_seed` is the magnitude used
    for the step. No rescale.
  - Under `step_size_policy = :adaptive_ffw`, `h_seed` is the *initial*
    magnitude fed to the FFW controller. The controller may then
    shrink it via `q ≤ 1` rescales but cannot grow it. The per-step
    accepted magnitude is therefore bounded above by `R(z_cur)`.
  - The walk termination distance — the predicate
    `abs(z_cur - target) ≤ ...` — uses `R(z_cur)` as well, so the
    walker approaches the target with the local-density resolution.

**Validation**: `R` is evaluated through a single helper
`_eval_node_separation(R, z, T)`. Non-finite or non-positive returns
throw `ArgumentError` with a suggestion. CLAUDE.md Rule 1 fail-loud.

## Alternatives considered

**A.  Wrapper grid-generator (option (a) from worklog 017).** A
higher-level function `nonuniform_grid(R, domain)` returning an
`AbstractVector{Complex}` consumed by the existing
`path_network_solve`. Rejected: splits ownership of the walker's
step-magnitude decision between the grid generator and the wedge
selector. FFW md:97 is explicit: the **step length** itself is
`|h(ζ)| = c·R(ζ)` — R controls the walker's spatial step, not just the
target-node placement. A wrapper that only shapes the targets and lets
the walker advance by uniform `h` between them recovers half the
benefit (correct target density) and forfeits the other half (correct
step density). The FFW algorithm is one indivisible thing.

**B.  Cell-list-driven Hjelle-Daehlen quasi-uniform sampling.** FFW
md:72 cite Hjelle & Daehlen 2006 as the source of their grid-generation
algorithm. Rejected for the reference impl: a Hjelle-Daehlen
implementation is ~500 LOC and replaces `R(ζ)` with a more capable but
more complex spatial-density specification (a callable that returns
both the target spacing AND a candidate-rejection criterion). For the
seven FFW figures in scope, the simple `R(ζ)` interface suffices; FFW
themselves use the same linear `R(ζ) = (8 - Re ζ)/20` for all PIII
figures and a similarly simple form for PV. File a separate bead if a
downstream user needs the full Hjelle-Daehlen capability.

**C.  Post-hoc tree pruning after a uniform Stage-1 walk.** Run the
walker with a small fixed `h` covering the densest pole-region
spacing, then prune the visited tree's smooth-region nodes to recover
a non-uniform distribution. Rejected: wastes the walker effort. The
densest-region `h` is `R(ζ = 8) = 0` in the FFW prescription; pruning
back from there means doing the most-expensive walk and throwing
~80% of it away.

**D.  Compose adaptive `h` only.** Skip R entirely and rely on
adaptive-FFW's `q ≤ 1` rescale to shrink the step in high-pole-density
regions. Rejected by experiment: adaptive shrinks `h` only when `T(h) >
Tol`, which is a **temporal** truncation-error signal, not a **spatial**
node-density signal. A smooth low-density region with a pole at
distance `5h` away has small `T(h)` until the walker is one step from
the pole — by which time the walker has already over-resolved a region
where there were no poles. R(ζ) provides the spatial-density target
*independent* of where the next pole happens to be; the two compose
orthogonally and the spec NU.1.7 pins this.

## Consequences

**Code**:
  - `src/PathNetwork.jl`: new `node_separation` kwarg (with `nothing`
    default) on `path_network_solve` and `_solve_with_schwarz_reflection`.
    New helper `_eval_node_separation(R, z, T)` validates R's output.
    Inline loop swap: `h_step` initial seed and walk-termination
    distance both use `R(z_cur)` when `node_separation` is supplied.
    Net effective LOC: 310 → 331 (+21). The file already exceeded the
    Rule 6 200-LOC cap pre-this-change (ADR-0011 ate ~70 lines for
    adaptive plumbing); the additional friction is documented in
    worklog 035 as a deferred split-point.
  - Module docstring extended with "Non-uniform Stage-1 node placement"
    section explaining the FFW prescription, the bounded-above
    composition with `:adaptive_ffw`, and the fail-loud validation
    contract.

**Tests**:
  - `test/non_uniform_nodes_test.jl` adds 7 testsets / 295 assertions
    (NU.1.1 - NU.1.7). All mutation-proven (M1-M4); footer logs bite
    counts (M1: NU.1.3 + NU.1.7 break; M2: NU.1.4 inverts; M3: NU.1.7
    bounded-above violations; M4: NU.1.3 + NU.1.7 break).
  - Aggregate count: 1751 → 2046 GREEN (+295).
  - NU.1.3 captures the headline: with the FFW Fig 1 prescription on
    a reduced ζ-plane patch (Re ζ ∈ [-1, 6], Im ζ ∈ [-2π/3, 2π/3]),
    the walker visits **517 nodes** (in-band: 80-3000). FFW report
    2701 for the full PIII Fig 1 domain (~4× wider strip + ~30% wider
    Re ζ range); 517 is consistent with the area ratio.

**Figures unblocked**:
  - FFW 2017 Figures 1-4 (PIII) and Figures 5-7 (PV) Stage-1 walks
    become reproducible. Composition with `SheetTracker` (Tier 5,
    already shipped) closes the full pipeline.

**Performance**:
  - Per-step cost under `node_separation` is one extra `R(z_cur)` call
    + one validation per inner loop iteration. The user-supplied `R`
    dominates if it's expensive; FFW's prescription is a single
    subtraction-and-divide, negligible.
  - Wall-time impact on a single FFW Fig 1 PIII patch (517 visited
    nodes, ~5000 inner step attempts): ~3s at order=30, Float64.
    Compare against ~30s for the equivalent uniform-`h = 0.05` walk
    that would resolve the same high-density region.

**Compatibility**:
  - Default `node_separation = nothing` preserves byte-identical
    output on every existing test, figure, and downstream caller
    (`Painleve`, `lattice_dispatch_solve`, etc.). Opt-in via kwarg.
  - The `enforce_real_axis_symmetry` Schwarz path threads
    `node_separation` through unchanged.

**Out of scope (deferred to separate beads)**:
  - `lattice_dispatch_solve` R(ζ) support — separate concern; file a
    bead when a downstream caller of the lattice dispatcher needs it.
  - BVP-layer non-uniform support — out of scope; BVP doesn't walk
    a 2D path-network.
  - Hjelle-Daehlen sampling — alternative B above; file a bead if a
    user needs richer spatial-density specifications than `R(ζ)`.

## References

  - FFW 2017 §2.1.2 — `references/markdown/FFW2017_painleve_riemann_surfaces_preprint/FFW2017_painleve_riemann_surfaces_preprint.md:67-72` (node separation rationale), md:97 (prescribed-step-size method), md:101 (Fig 1 2701-point pin).
  - FW 2011 §3.1 — `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:155-166` (the FW 2011 wedge walker that this enhances).
  - ADR-0004 — path-network architecture (consumed).
  - ADR-0011 — adaptive Padé step (composed orthogonally; see "Composition" above).
  - `docs/worklog/017-coord-transforms-pIII-pV.md` §"What is NOT shipped" — the prior deferral note.
  - `docs/worklog/035-non-uniform-stage-1-nodes.md` — implementation diary.
  - `src/PathNetwork.jl` — module docstring section "Non-uniform Stage-1 node placement (`node_separation`, opt-in)" + the inline R-driven seed in `path_network_solve`.
  - `test/non_uniform_nodes_test.jl` NU.1.1-NU.1.7 — acceptance tests + mutation-proof footer.
