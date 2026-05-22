# ADR-0026 — Vector path-network: a resilient walk + dense targets for the headline figure

**Status**: Accepted (2026-05-22). Supersedes the mis-framed deferred
bead `padetaylor-0ln.38` ("2D re-expansion lattice — a different
architecture"). **Bead**: `padetaylor-0ln.40`. **Follows**: worklog 059
(the honest reassessment), ADR-0025 Amendment 13.

## Context

Worklog 059 correctly judged the P_I⁽²⁾ tritronquée headline figure
inadequate: its pole-rich wedge is ~95 % blank. It then diagnosed the
cause as "`src/VectorPathNetwork.jl` is a minimal ~156-LOC port of
v0.1's 1108-LOC `src/PathNetwork.jl`; the region-filling machinery —
resilient walk, fine-grid targets, **barycentric Stage-2 fill** — was
never ported," and framed the fix (bead `0ln.40`) as "port the full FW
path-network driver."

A four-agent deep-dive recon (2026-05-22; the `0ln.37` plan-first
pattern) read every relevant file and the FW 2011 markdown end-to-end.
**Two of worklog 059's four diagnosed gaps do not exist.** This ADR
records the corrected diagnosis and scopes the actual work.

## The corrected diagnosis

| Worklog 059 claim | Recon finding |
|---|---|
| Missing the **barycentric Stage-2 fill** | There is no barycentric fill to port. The scalar `PathNetwork.jl:669–693` Stage-2 is *also* plain nearest-node single-Padé lookup. FW 2011 itself prescribes exactly that — "one single Padé step suffices … no tests needed" (FW-md:166, :395–397); the word *barycentric* appears in FW only for BVP Chebyshev post-processing. A multi-node distance-weighted blend *was* built and auditioned on the vector side (the "B4 audition") and **rejected on measured evidence** (`VectorPathNetworkStage2.jl:133–183`). Nothing to port; nothing missing. |
| Missing the **Stage-2 fine-grid fill** | It was ported — `VectorPathNetworkStage2.jl` exists (bead `0ln.20`, commit `1d3c32b`), wired through the `fine_grid` kwarg. Worklog 059 is stale here. |
| **Brittle walk** — `throw`-aborts on first failure | **True** — and the scalar `PathNetwork.jl` is *equally* brittle (it has no skip-and-continue either). Resilience is genuinely *new* work for both drivers, not a port. |
| **171 coarse targets**, a fill needs ~25 000 | **True.** The solver already accepts any target vector verbatim (`VectorPathNetwork.jl:428–442`, `:498`); the figure helper simply passes a sparse 9×19 fan. A denser target set is a *caller-side* change. |

The vector path-network is therefore **not a skeleton**. It already
carries machinery the scalar driver lacks — the B1 true-radius Stage-2
gate and the B2 `:max_q_root` adaptive wedge walk (ADR-0025). A
wholesale 1108-LOC port would *duplicate and in places regress* it.

The figure is blank for exactly two reasons: **(1)** the Stage-1 walk
`throw`s and aborts the entire run on the first unreachable target, and
**(2)** the figure asks for only 171 filament-endpoint targets, so the
path-network tree is 171 filaments with wide gaps between them rather
than a space-filling tiling.

## Decision

Three binding decisions. No 1108-LOC port; no Schwarz reflection (the
user has ruled it out, and FW 2011 uses conjugate symmetry only as an
accuracy *diagnostic*, FW-md:303–310, never as a work-halving shortcut);
no barycentric Stage-2 (it does not exist in FW and was already
evidence-rejected here).

### D1 — A resilient Stage-1 walk (new capability, `src/VectorPathNetwork.jl`)

`vector_path_network_solve` gains an additive kwarg
`on_target_failure::Symbol = :throw`:

- `:throw` (default) — byte-identical to today. Every existing call
  site and the ~2600 existing assertions are unaffected.
- `:skip` — the per-target walk is wrapped so that a failure (any of
  the three throw sites: unreachable `VectorPathNetwork.jl:515`,
  all-five-candidates-fail `VectorWedgeStep.jl:308`, adaptive-step
  collapse `VectorWedgeStep.jl:380`) **does not abort the run**. The
  failing target is recorded and the walk continues to the next target.

This is **not** a silent swallow — that would violate Rule 1. It is
fail-*loud-by-accounting*: every skip is captured as first-class data.
`VectorPathNetworkSolution` gains a tenth field `failed_targets`, a
concretely-typed vector of records, one per skipped target, carrying:
the target point, a `reason::Symbol`
(`:unreachable` | `:all_candidates_failed` | `:step_collapse`), the
`z` at which the walk stuck, and the residual distance `|z − target|`.
Backward-compat constructors (6-/8-/9-arg) default it to empty, exactly
as `visited_jets` was added (`VectorPathNetwork.jl:298–315`).

Partial nodes laid down by a walk that later failed are **kept** — they
are valid solution nodes (the walk pushes to `visited_*` only on an
*accepted* step) and still contribute coverage. A failed walk is "this
target was not reached," never "discard the filament."

The contract: with `on_target_failure = :skip`, the solver returns
normally; a non-empty `failed_targets` is a signal the *caller must
surface* (the figure reports the count and the achieved coverage).

### D2 — Dense, disc-spaced wedge targets (caller-side, the figure helper)

The figure's wedge target generator (`surf_wedge_targets` /
`SURF_TARGET_RADII`/`SURF_TARGET_ANGLES` in
`figures/_kkg_pi2_surface_helpers.jl`) is replaced by a generator that
**tiles the `|arg x| < 36°` wedge with targets at disc spacing**, so
adjacent path-network filaments' validity discs overlap. FW's analogue
is its 40×40 = 1600-node coarse grid feeding a 161×161 fine grid
(FW-md:153–164). The B1 true-radius gate measures honest disc radii at
median ≈ 0.06, p90 ≈ 0.16 (ADR-0025); target spacing is therefore a
*measured* parameter — start at ≈ 0.25 (≈ 2× the disc radius), render,
**measure the achieved coverage**, and tighten until coverage saturates
or runtime forces a stop. The expected envelope is a few thousand
targets; the realised count and coverage are reported, not assumed.

### D3 — Walk failures are investigated through the FW lens

A skip is a safety net, not an answer. Per CLAUDE.md Rule 2 ("all bugs
are deep"), every distinct failure cluster surfaced by `failed_targets`
is **root-caused before the figure is called done**, with the explicit
question: *how did FW solve the analogous problem?* FW's random-path
strategy empirically reached "within h of all target points" (FW-md:163)
— FW reports essentially **zero** unreached targets. So a vector walk
that skips a non-trivial fraction is exhibiting a *bug*, and FW's
near-100 % success is the benchmark. Candidate FW-side analogues to
check against: the §3.1 nearest-visited-node restart, the §5.5
"complete the pole field before stepping into smooth regions" path
ordering, the (order, h) = (30, 0.5) standard. If, after investigation,
a residual frontier is genuinely intractable, it is **measured and
reported honestly** (the ADR-0025 Amendment-13 "numerically-supported
frontier" doctrine), never quietly deferred.

## Out of scope

- The wholesale 1108-LOC port of `PathNetwork.jl` (recon shows it is
  unnecessary and partly regressive).
- Schwarz / conjugate-reflection symmetry (user-ruled-out; FW-faithful
  to omit).
- A barycentric / multi-node Stage-2 blend (absent from FW; already
  evidence-rejected by the B4 audition).
- Non-uniform `node_separation`, branch-cut tracking, the diagnostics
  layer — the other vector-port deferrals; none is needed for a filled
  honest wedge.
- The pole-free sector (the triple-method majority vote, ADR-0024) —
  sound and unchanged.

## Verification / acceptance

1. The B2 `:max_q_root` adaptive walk, the B1 true-radius Stage-2 gate,
   and the VC-4…VC-10 validation suite (`figures/test_kkg_pi2_surface.jl`)
   are **unchanged and stay GREEN**.
2. New resilience behaviour is TDD'd in `test/vector_path_network_test.jl`
   and **mutation-proven** (perturb the skip/accounting logic, confirm
   RED, restore).
3. `on_target_failure = :throw` default ⇒ the full `Pkg.test()` suite
   stays GREEN with no assertion changes (additive kwarg + additive
   struct field).
4. The headline figure re-renders with a **genuinely filled** pole
   wedge; the achieved honest coverage fraction is **measured and
   reported** in worklog 060 and in the figure's docstring.
5. Any residual unreached frontier is investigated per D3 and reported.

## Child beads

- `0ln.40.a` — resilient walk: `on_target_failure` kwarg +
  `failed_targets` field + accounting (D1). TDD + mutation-proof.
- `0ln.40.b` — dense disc-spaced wedge target generator (D2).
- `0ln.40.c` — re-render the headline figure; measure & report coverage.
- `0ln.40.d` — FW-lens investigation of walk failures + remediation (D3);
  conditional on what `.c` surfaces.
- `0ln.40.e` — re-verify (VC suite + full `Pkg.test()`), worklog 060,
  ADR finalisation, epic bookkeeping.

## References

- `docs/worklog/059-headline-figure-honest-reassessment.md`
- `docs/adr/0025-headline-figure-re-resolution.md` (Amendment 13)
- `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md`
  — §3.1 path network (:151–168), §5.4.2 Stage 2 (:395–397),
  conjugate-symmetry diagnostic (:303–310)
- `src/VectorPathNetwork.jl`, `src/VectorWedgeStep.jl`,
  `src/VectorPathNetworkStage2.jl`, `src/PathNetwork.jl`
