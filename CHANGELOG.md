# Changelog

All notable changes to PadeTaylor.jl are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
the project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

Post-v0.1.0 work, all on `main` and pushed.  Test suite **9369 pass / 2
broken / 0 fail**, re-confirmed by `scripts/quality_gate.sh full` on 2026-08-23
(Julia 1.12.3, 17m50s; unchanged since the 2026-06-30 Fig 4.7 seam cure; up from
1311 at the v0.1.0 tag).  44 source modules
(up from 14); 34 ADRs (up from 4); 79 worklog shards (up from 18).

> The two `@test_broken` markers are intentional deferred-feature fixtures
> (`cm-n2` in `corpus_vector_polefield_test`, `pi2-tritronquee` in
> `corpus_painleve_rational_test`) that auto-flip to "Unexpected Pass" the
> day their bug is fixed.  Investigate only FAILs.  See
> `scripts/quality_gate.sh` "EXPECTED-NOISE".

### Tests — BigFloat-256 tritronquée pin (bead `padetaylor-7zw`)

- **`padetaylor-7zw`** — `test/fw_fig_41_test.jl` gains a parallel
  BigFloat-256 testset (8 assertions) rerunning the FW Fig 4.1 step-(i)
  BVP at N=240 with the FW eq. 4.1 targets as `big"…"` literals
  (FW2011 md:226).  Pins rel-err u(0) ≤ 1e-11 (measured 1.85e-12),
  rel-err u'(0) ≤ 3e-10 (measured 1.75e-10), and the 256-bit floors
  ‖R‖_∞, |imag| ≤ 1e-60 (measured 1e-68 / 1e-76) as the precision
  witness.  Finding: at N=240 the u(0)/u'(0) error equals the Float64
  figure — truncation-limited, not precision-limited; N=320 reaches FW's
  16-digit limit (1.6e-16) in 75 s.  Mutation-proved (M1 reference
  digit, M2 Float64 leak).  Worklog 080.  ≈60 s added to the suite.

### Fixed — quality gate

- `scripts/quality_gate.sh fast` was RED on a clean tree: it ran
  `test/quality_test.jl` under `--project=.` where the `[extras]` test deps
  (Aqua, ExplicitImports) do not resolve; it had only ever passed via the
  user's global env. Now `Pkg.test(test_args=["quality-only"])` with a
  short-circuit in `runtests.jl` (bead `padetaylor-0dw7`).

### Fixed / Tests — validation audit 2026-08-23 (beads `padetaylor-lkrk`, `padetaylor-cwcp`, `padetaylor-044m`)

- **`padetaylor-lkrk`** — five input-validation gaps closed, each RED-first
  with a Rule-1 `ArgumentError` + suggestion and mutation-proven (guard
  removed → RED, restored): (1) `edge_gated_pole_field_solve` now checks
  EVERY axis step (`all(diff(xs) .≈ Δx)`, same for `ys`), as its docstring
  always promised, not just the first; (2) `max_iter ≥ 1` is required and
  exhausting it before the region-growing fixpoint THROWS with the field /
  frontier sizes and pass count instead of returning a half-grown field;
  (3) `step_jorba_zou` / `vector_step_jorba_zou` require jet length ≥ 3
  (`p ≥ 2`) — a length-2 jet put `k = 0` in the candidate set and
  `(ε/|c₀|)^(1/0)` silently yielded `h = 0`; (4) `edge_gated_pole_field_solve`
  and `windowed_path_network_solve` promote integer / mixed-precision axes
  up front to `float(promote_type(eltype(xs), eltype(ys)))` — RED-first
  showed a `Float32` `xs` with `Float64` `ys` silently DEMOTED the whole
  solve to single precision (the audit's predicted `TypeError` for `Int`
  axes did not reproduce: `new` converts ranges); (5) `EdgeGatedSolution.
  field_mask` docstring now says "admitted/anchored region incl. the
  unclassified seed disc", not "confirmed by the edge detector".  Tests:
  EG.2.1 (+11), EG.2.2 (new), 4.1.6 (new), VSC.1.5 (+2), FSEAM.4 (new).
- **`padetaylor-cwcp`** — `src/LinAlg.jl` "Why we expose `full`" docstring
  said the null right singular vector is `transpose(Vt[end,:])`; `Vt = V'`
  is the conjugate transpose, so for complex input it is `conj(Vt[end,:])`
  (code unaffected — `RobustPade.jl:450` takes `abs.()` and recomputes `b`
  by adjoint-QR).  New `linalg_test.jl` 1.1.3c: complex rank-deficient 3×4
  with a genuinely complex null vector, both backends (LAPACK `ComplexF64`,
  GenericLinearAlgebra `Complex{BigFloat}`); asserts `‖A·conj(Vt[end,:])‖ ≤
  1e-12‖A‖` AND that the unconjugated row is NOT null (`> 0.1‖A‖`).
  Mutation `return conj(F.Vt)` → 6 RED in 1.1.3c, 1.1.3b (real) GREEN.
- **`padetaylor-044m`** — `test/robustpade_test.jl` header over-claimed
  "six captured oracles"; now lists per case exactly what is pinned (2.1.1 /
  2.1.2 / 2.1.5 coefficients; 2.1.3 / 2.1.6 degrees only; 2.1.4 degrees +
  poles).  2.1.4 now asserts the captured `padeapprox.m` tan(z⁴) pole set
  (`_oracles.jl:28-31`, `padeapprox.m:150`) as a bidirectional set match —
  measured max nearest-neighbour distance `1.35e-14`, pinned `1e-12`; the
  inner ring is anchored to the true `(π/2)^(1/4)` to `1e-5` (measured
  2.1e-6), the outer ring is an approximant artefact (1.4993 ≠ `(3π/2)^(1/4)`
  = 1.4734) and is pinned to the capture only.  The 2.1.7 note wrongly said
  2.1.3 carries a coefficient match — corrected.  New mutant in
  `test/mutation/run_mutation_gate.jl` (catalogue #6, 13 total): drop the
  `padeapprox.m:111-117` QR-reweighting (`src/RobustPade.jl:450-462`) —
  hand-verified 5 RED (2.1.2 coefs ×2, 2.1.4 ν + pole set ×2), gate run on
  index 6 alone → KILLED.

### Tests — rigour audit 2026-08-23 (beads `padetaylor-aqw1`, `padetaylor-urmg`, `padetaylor-ll3t`)

- **`padetaylor-aqw1`** — `test/pathnetwork_test.jl` PN.2.4 pins FW 2011
  Table 5.1 column (c), `u(28.261) = 9.876953517025014e6` (end point high on
  a pole wall; `references/markdown/FW2011_painleve_methodology_JCP230/
  FW2011_painleve_methodology_JCP230.md:370-372`, table at `:385-391`).
  The oracle constant had sat in `test/_oracle_problems.jl` since Phase 6
  with no assertion.  Measured F64 rel-err `6.88e-10` (FW's own Padé:
  `7.92e-10`); pinned `rtol = 5e-9` plus a relative imag bound.  BF-256
  one-shot `7.83e-11` (not in-suite, ~52 s).  Mutation-proven: 8th-digit
  perturbation of the oracle → RED.
- **`padetaylor-urmg`** — the harmonic `[cos, −sin]` pins in
  `test/vector_step_control_test.jl` (VSC.1.4) and `test/vector_problems_test.jl`
  (VP.1.2, VP.1.4) tighten `atol 1e-8 → 5e-13`.  Re-measured max error
  `5.1e-14` (was `1.27e-9` at `decc32d`); bisected to commit `8f2cbfd`, the
  ADR-0028 dispatch wiring, exactly as the now-deleted "forthcoming ADR-0028
  will tighten" comment predicted.  Mutation-proven: forcing cell A in
  `shared_pade_select` → 21/69 and 26/72 RED.
- **`padetaylor-ll3t`** — `test/heun_test.jl` gains the MUTATION-PROOF footer
  it lacked: sign-flipping the lagged term of the HeunG (`_heun_g_frobenius_at_0`)
  and HeunC (`_heun_c_frobenius_at_0`) recurrences each RED the corresponding
  30-digit oracle testset (HG.2 3/6, HC.2 3/6); restored, 49/49 GREEN.

### Tests + docs — rigour audit 2026-08-23, second batch (beads `padetaylor-ps4g`, `padetaylor-divo`, `padetaylor-fqwf`)

- **`padetaylor-ps4g`** — `test/extrapolate_test.jl` SX.1.2b value-pins
  `eval_at(; extrapolate=true)` against the entire closed form `exp(z - z₀)`
  at six points outside every visited disc (local `|t| ∈ (1, 3]`; default mode
  asserted NaN there).  Measured rel. err ≤ 7.2e-16 for `u` and `u'`; pinned
  at `rtol = 1e-12`.  Far-out degradation recorded (7.0e-6 at `z = 10`).
  Mutation X4 (wrong local `t`): 12 RED of 48; the old `isfinite` checks alone
  did not bite.
- **`padetaylor-divo`** — `src/VectorStepControl.jl` docstring and ADR-0021
  corrected: the "2-norm step is never larger than the ∞-norm step" ordering
  holds only in absolute mode; in relative mode `ε` is rescaled by the same
  `vnorm(c₀)` and the 2-norm step can exceed the ∞-norm step by up to `√d`
  (Jorba–Zou eq. 10 uses the sup norm, md:601-604).  Policy stays 2-norm.
  New test VSC.1.3b pins the counterexample ratio `h₂/h_∞ = √2`
  (measured 1.4142135623730954); mutation M4 (ε via ∞-norm regardless of
  `vnorm`): 4 RED of 73.
- **`padetaylor-fqwf`** — `windowed_path_network_solve` now `@warn`s when
  `window_extent ≥` an axis's span (1×1 tile grid on that axis: the composite
  silently degenerates to the monolithic path-dependent solve), naming the bead
  and the fix (reduce `window_extent`); documented in the docstring validation
  list; ADR-0034 no longer claims a warning that did not exist.  New test
  FSEAM.3 (`@test_logs`, both warn and must-be-silent cases); mutation F3
  (never warn): 2 RED of 7.  The warning pushed `src/WindowedComposite.jl`
  to 211 effective LOC, so the tiling geometry (`_tile_centers`,
  `_nearest_center`, `_window_seed`, `_warn_single_window`) moved verbatim
  into new `src/WindowedTiling.jl` (50 LOC, literate docstring citing
  ADR-0034 + FW 2011 md:141-147); `WindowedComposite.jl` is now 166 LOC.
  Pure relocation — no behaviour change.

### Fixed — input validation in the Padé layer (beads `padetaylor-ked0`, `padetaylor-lbqb`)

  - **`src/SharedPade.jl`** — the zero-input guard was absolute (`‖c‖₂ ≤ tol`) and falsely rejected valid small-amplitude jets; it is now `iszero(‖c‖)`, restoring the scale covariance every other threshold already had (GGT 2013 `τ = tol·‖c‖₂`, md:213/225). Pinned by `SP.6` (α ∈ {1e-20, 1e20}); mutation-proven (`padetaylor-ked0`).
  - **`src/RobustPade.jl` + `src/SharedPade.jl`** — `robust_pade` / `shared_denominator_pade` now throw `ArgumentError` ("the Taylor jet overflowed; reduce h or order") on any `Inf`/`NaN` coefficient or `tol ∉ [0, Inf)` instead of silently returning the zero approximant. Pinned by `2.1.9` and `SP.7`; mutation-proven (`padetaylor-lbqb`).

### Fixed — Fig 4.7 / 4.8 path-dependence seam (bug `padetaylor-vwgl`, ADR-0034, worklogs 077/078/079)

  - **`src/WindowedComposite.jl` — bounded-window composite solve** cures the
    seed-dependent grain boundary in `path_network_solve`'s monolithic walk tree.
    Tiles the domain into overlapping windows, solves each independently from the
    shared IC at a per-window seed (`_window_seed(rng_seed, wi)`), and composites
    by Voronoi-core ownership.  Confirmed numbers (5×5 composite over `[-50,50]²`,
    two global seeds): pole-count seed-variance |Δ| **151→12**, pole-set match
    **77.1%→99.4%**, median pole displacement **0.076→0.000**, ~2× faster than
    monolithic.  New public API: `windowed_path_network_solve`,
    `windowed_extract_poles`, `WindowedCompositeSolution`.
  - **`edge_gated_windowed_poles`** — production driver for Fig 4.7 / 4.8: run
    `edge_gated_pole_field_solve` for the bloom-free pole-field mask, run the
    windowed composite, filter extracted poles to the mask.  Measured on the PI
    tritronquée over `[-50,50]²`: off-wedge bloom **4475→0**; 97.6% agreement
    with the edge-gated baseline in the pole field.  `fw2011_fig_4_7.jl` and
    `fw2011_fig_4_8.jl` rewired to this driver.
  - **`test/field_seam_test.jl`** (bead `padetaylor-sny7`) — non-gameable
    seed-invariance gate: FSEAM.1 asserts two global seeds genuinely re-randomise
    every window's walk; FSEAM.2 asserts the composited pole field is
    seed-invariant (bidirectional match ≥ 0.90 — composite 97.0–97.2% GREEN,
    monolithic 77.6% RED).  Mutation-proven: M-frozen-seed → FSEAM.1 RED,
    M-monolithic → FSEAM.2 RED (both executed + restored byte-clean).
  - **INGN ℘ verification** (bead `padetaylor-ingn`, probe
    `p6_ingn_pole_truth.jl`, Weierstrass-℘ oracle at `HALF=30`): windowed
    recall=1.000 == monolithic (drops zero real poles); windowed
    precision=1.000 > monolithic 0.998 (suppresses spurious seed-dependent
    poles).  The composite's ~220 fewer poles are spurious ones correctly
    removed, not real poles dropped.
  - **Residual** (Rule 9, bead `padetaylor-us19`, P3): the interior pole lattice
    (≥3 rings inside the mask) is seed-invariant at 97.6–97.8% (seam cured);
    perimeter poles flicker at 79.9–89.6% match (boundary-pole flicker, not
    blocking, benign for single-seed figures).

### Fixed — out-of-class guard (bug `padetaylor-v1ub`, ADR-0033)

  - **`solve_pade` now FAILS LOUD on out-of-class (non-meromorphic)
    input** instead of silently returning finite, plausible, *wrong*
    values (the only confirmed silent-wrong-answer bug).  On a recast
    whose solution has an essential singularity (`u'' = u(1+2z)/z⁴`,
    exact `u = e^{1/z}`) the solver used to integrate toward `z = 0`
    returning finite values that degrade to relerr ≈ 1.2e17 with the
    WRONG SIGN, no throw, no NaN — a Rule-1 violation.  The new guard
    watches the per-step **two-order Padé convergence defect** δ (the
    disagreement between the working `[m,n]` and reduced `[m-1,n-1]`
    Padé of the same scaled jet at the step endpoint) and throws an
    `OutOfClassError` once δ exceeds τ = 1e-3 AND has grown monotonically
    over the last K = 2 steps — the de Montessus / Nuttall–Pommerenke
    signature of a jet leaving the meromorphic class (GGT 2013 §8).  The
    history gate makes the guard safe: it fires on the e^{1/z} approach
    but CANNOT trip on the single-step across-0 bridge or on legitimate
    pole bridging (δ stays at the rational-approximation floor while
    bridging a pole — measured max in-class δ ≈ 4e-7, ~3.8 orders below
    τ).  Default-on; pass `check_in_class = false` to opt out (legacy
    behaviour, byte-for-byte identical numerics).  New module
    `src/OutOfClass.jl` (literate); the perf-critical
    `PadeStepper.pade_step_with_pade!` hot path is left byte-for-byte
    unchanged.  See ADR-0033 and `test/corpus_out_of_class_test.jl`.

### Added — worklogs 020–033

  - `RobustPade.classical_pade_diagonal` — the classical FW 2011
    §5.1.4 Toeplitz-backslash diagonal Padé.  Now the default at
    `Float32` / `Float64` / their `Complex` variants (faster and more
    accurate on smooth inputs than the SVD path); the GGT 2013 SVD
    route remains the default for `BigFloat` / `Arb`.  See ADR-0005.
  - `PoleField.extract_poles` — reads pole locations back out of a
    solved `PathNetworkSolution` / `PadeTaylorSolution` by rooting the
    stored Padé denominators and clustering across nodes.
  - `EdgeGatedSolve.edge_gated_pole_field_solve` + `EdgeGatedSolution`
    — region-growing pole-field solver that confines the IVP walk to
    the pole field (morphological open + flood-fill), curing the
    spurious-pole bloom of plain `path_network_solve` on solutions
    with large smooth sectors (FW 2011 line 401).
  - `Painleve.PainleveProblem` (`src/Painleve.jl`) — per-equation
    problem builder for all six Painlevé equations (PI/PII/PIV built
    directly; PIII/PV/PVI via the existing coordinate transforms).
    See ADR-0006.
  - `Painleve.PainleveSolution` (`src/PainleveSolution.jl`) —
    self-describing solve-output wrapper carrying the Painlevé
    identity + coordinate frame, with a uniform `z`-frame access
    surface (`sol(z)`, `poles`, `grid_values`, `equation`,
    `parameters`, `solutionname`).  Makes the `:transformed`
    equations (PIII/PV/PVI) work end-to-end through
    `path_network_solve`.  See ADR-0007.
  - `tritronquee` / `hastings_mcleod` (`src/PainleveNamed.jl`) —
    named-transcendent constructors with literature-pinned 16-digit
    initial conditions baked in (PI tritronquée from FW 2011 §4.1;
    PII Hastings–McLeod from FW 2014).  See ADR-0008.
  - `PadeTaylorMakieExt` — Makie plot recipe (`painleveplot`) for
    `PainleveSolution`.
  - `figures/` — standalone Julia project reproducing **thirteen**
    FW 2011 figures as runnable scripts (Fig 3.1–3.3, 4.1, 4.2–4.4,
    4.7, 4.8, 5.1, 5.2), each writing a PNG to `figures/output/`.
  - Documenter.jl docs site at `docs/build/`.
  - ADRs 0005–0010: classical-Padé default at `Float64` (0005),
    `PainleveProblem` layer (0006), `PainleveSolution` wrapper (0007),
    named-transcendent constructors (0008), `EdgeDetector` h-aware
    level (0009), `pii_rational` / `pii_airy` / `piv_entire`
    closed-form Painlevé families (0010).

### Added — FFW 2017 arc (worklogs 034–045)

The 11-step plan to reproduce the seven FFW 2017 figures.  **8 of 11
steps shipped** as of worklog 045.

  - **A1 adaptive Padé step** (`:adaptive_ffw` step_size_policy on
    `path_network_solve`, worklog 034, ADR-0011) — FFW 2017 §2.1.2
    truncation-error controller `q = (k·Tol/T(h))^(1/(n+1))`.
  - **A2 non-uniform Stage-1 nodes** (`node_separation::Function`
    kwarg, worklog 035, ADR-0012) — FFW Fig 1 prescription
    `R(ζ) = (8 - Re ζ)/20`.
  - **A3 η-plane PVI** (`pVI_eta_transformed_rhs` + `pVI_z_to_η` /
    `pVI_η_to_z` + `:transformed_eta` frame on `PainleveProblem(:VI)`,
    worklog 041) — FFW 2017 eq. 5 (md:154) verbatim; branch-point-free
    region `Re η < log(2π)`.
  - **A4 constrained-wedge routing + per-branch sheet bookkeeping**
    (new `src/BranchTracker.jl` module with `segment_crosses_cut`,
    `any_cut_crossed`, `step_sheet_update`, `resolve_cut_angles`;
    `branch_points` + `branch_cut_angles` + `cross_branch` +
    `initial_sheet` kwargs on `path_network_solve` + new
    `visited_sheet` field on `PathNetworkSolution`, worklog 042,
    ADR-0013) — FFW 2017 §2.2.2 (md:163-189).
  - **A5 sheet-aware Stage-2** (`grid_sheet` kwarg on
    `path_network_solve` + new public `eval_at_sheet(sol, z, sheet)`
    accessor, worklog 043) — Stage-2 lookup restricted to matching-
    sheet visited nodes.
  - **A6 IVP+BVP hybrid driver** (`solve_pole_free_hybrid` +
    `IVPBVPSolution` + `pIII_asymptotic_ic` in new
    `src/IVPBVPHybrid.jl`, worklog 039, ADR-0014) — FFW 2017 §3
    PFS↔BVP coupling for pole-free sectors.  Additive 3-arg-RHS
    overload `bvp_solve(f, ∂f_∂u, ∂f_∂up, ...)` for the `(w')²/w` PIII
    term.
  - **`extrapolate=true` Stage-2 kwarg + new public `eval_at(sol, z;
    extrapolate=false)` accessor** (worklog 045, ADR-0015) — aligns
    Stage-2 with FFW md:62 spec (evaluate Padé at every fine-grid
    cell regardless of disc radius).  Default `false` preserves
    CLAUDE.md Rule 1 fail-soft NaN contract; opt-in `true` fills
    figure renders without white gaps.
  - **B1 FFW Fig 6** PV generic three-sheet (`figures/ffw2017_fig_6.jl`,
    worklog 036) — first FFW 2017 figure reproduced; per-sheet errors
    beat FFW by 2-3 orders.
  - **B2 FFW Fig 1** PIII three-sheet spiral (`figures/ffw2017_fig_1.jl`,
    worklog 037) — FFW's headline figure; sheet-0 conjugate-symmetry
    median 4e-15 beats FFW Exp-2's 1e-6 by 9 orders.
  - **B3 FFW Fig 4** PV tronquée three-sheet
    (`figures/ffw2017_fig_4.jl`, worklog 038).
  - **B4 FFW Fig 5** PIII tronquée + cond-number heatmap
    (`figures/ffw2017_fig_5.jl`, worklog 040) — uses A6 hybrid driver;
    cond-number pin `κ_r(z=30) ≈ 157` matches FFW md:264.
  - **B5 FFW Fig 2** PVI three-method reproduction (η + ζ-refuse +
    ζ-cross, `figures/ffw2017_fig_2.jl`, worklog 044) — first
    end-to-end demo of A1+A2+A3+A4+A5+A6 stack composition.

### Changed

  - All four updated FFW 2017 figure scripts (Fig 1, 2, 4, 6) opt into
    `extrapolate=true` Stage-2 and render at **100% coverage** (was
    0.9–98% pre-ADR-0015).

### Added — Diagnostics layer (worklog 048, bead `padetaylor-5t4`)

  - `Diagnostics.quality_diagnose(sol) → DiagnosticReport` — first-class
    loop-closure quality certificate on `PathNetworkSolution`.  Computes
    the per-edge midpoint disagreement `ΔP_rel := |P_A(M) - P_B(M)| /
    (|P_A(M)| + |P_B(M)| + ε)` over every non-tree Delaunay edge on
    sheet 0 of the visited-node graph, returning a `DiagnosticReport`
    with `(n_nodes, n_nontree, median, p90, p99, max, n_above_tol,
    n_catastrophic, edge_reports)` fields.  Motivated by the trimodal
    distribution found in the Fig 1 loop-closure probe
    (`external/probes/loop-closure-fig1/REPORT.md:34-77`): 6.3 % of
    loop closures catastrophic (`ΔP_rel > 1e-3`), clustering at
    high-Re ζ and near sheet boundaries.  See ADR-0016.
  - `PathNetworkSolution` gains `diagnostics::Union{Nothing,
    DiagnosticReport}` field (11th); backward-compat constructors
    preserve all existing call sites.
  - `path_network_solve` gains `diagnose::Bool = false` kwarg;
    `diagnose=true` attaches the report to the returned solution.
  - `PadeTaylorDiagnosticsExt` — weak-dep extension loading
    `DelaunayTriangulation` (ADR-0016; same pattern as ADR-0003).
    Without the extension, `quality_diagnose` throws a helpful
    `ErrorException`.
  - +32 new DG.* assertions across 8 testsets in
    `test/diagnose_test.jl`; mutation-prove A applied (RED→GREEN
    verified, 18 s, single-file isolation — full `Pkg.test()` deferred
    due to OOM friction, see worklog 048).  Existing 2508 assertions
    unchanged.

### Fixed — Lattice-dispatcher v3 (worklog 049, bead `padetaylor-0tj`)

  - **`lattice_dispatch_solve` BC-corruption divergence** — on PI tritronquée
    at large N (e.g. 121×121 over `[-20,20]²`, `h_grid = 0.333`), the plain
    `path_network_solve` IVP source walked into smooth pole-free sectors and
    returned corrupted values there (FW 2011
    `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:401`:
    "smooth regions are unstable regions for any IVP solver").  Those corrupted
    values became the Dirichlet BCs for the per-row BVP, causing Newton
    divergence (`‖Δu‖_∞ = 0.4458` vs `tol = 1.8e-12`).  Fixed by switching
    the default IVP source to `edge_gated_pole_field_solve` (ADR-0017;
    bead `padetaylor-0tj`).  Smooth cells outside the gated field are
    `NaN + NaN·im` — Rule 1 honesty.  Manual `mask` kwarg still routes through
    plain `path_network_solve` (API contract preserved).

### Changed — Lattice-dispatcher v3 (worklog 049, bead `padetaylor-0tj`)

  - **`lattice_dispatch_solve` default IVP source** — switched from plain
    `path_network_solve` to `edge_gated_pole_field_solve` (ADR-0017;
    bead `padetaylor-0tj`).  Consequence: smooth cells outside the gated field
    are `NaN + NaN·im` by default (was: plausible-but-wrong IVP-integrated
    values).  Manual `mask` kwarg still routes through plain `path_network_solve`.
  - **LD.1.1 pre-existing assertion adjusted** — the LD.1.1 `bvp_solutions == 0`
    assertion (which expected zero BVP cells under the old full-grid IVP default)
    is updated to `bvp_solutions ≤ 5`; under the new edge-gated default the
    one-ring dilation provides a small finite number of valid BVP flanks.
    The semantic flip is documented in `test/lattice_dispatcher_test.jl`.

### Added — Lattice-dispatcher v3 (worklog 049, bead `padetaylor-0tj`)

  - **`strict::Bool = true` kwarg on `lattice_dispatch_solve`** — `strict = false`
    catches Newton non-convergence (`ErrorException` matching the exact message
    from `src/BVP.jl:332-338`) and tags the affected cells `:bvp_fail`
    instead of throwing.  `strict = true` (default) preserves the v1/v2
    fail-fast behaviour; every existing test invariant is unchanged.
  - **`:bvp_fail` added to the `region_tag` enum** — was `{:ivp, :bvp,
    :ivp_only}`; now `{:ivp, :bvp, :bvp_fail, :ivp_only}`.  Downstream
    consumers pattern-matching on the enum need a fourth arm (none currently do).
  - +61 new LD.X.* assertions across 7 testsets in
    `test/lattice_dispatcher_test.jl`; mutation-X applied (RED→GREEN verified,
    23.4 s, single-file isolation — full `Pkg.test()` deferred due to OOM
    friction, see worklog 049).

### Added — API audit (worklog 050, bead `padetaylor-qur`)

  - **`docs/api-review-2026-05-16.md`** (743 LOC) — holistic pre-v1.0
    review of the 39 exported symbols across 23 modules: kwarg naming,
    return-type shape, error-message style, `!`-suffix discipline,
    PascalCase/snake_case consistency, export-list completeness.  Advisory
    only; no source files modified.  Identified 3 v1.0-blocking findings
    (canonical `h` rename → `padetaylor-xds`; `max_iter` normalisation →
    `padetaylor-0xn`; `pii_*` / `piv_entire` export omission →
    `padetaylor-gvz`) + 5 P3 cosmetic + 4 P4 defer.  12 follow-up beads
    spawned in total; these constitute the v1.0 normalisation work backlog.

### Changed — public-API kwarg canonicalisation (beads `xds` + `0xn`)

  Implements the three v1.0-blocking findings of the API audit
  (`docs/api-review-2026-05-16.md` §3(a).1, §3(a).2, §9 rows 1-2).
  All renames are **backward-compatible** via deprecation shims (see
  the "Deprecated" section below); every in-repo call site (src, tests,
  figures, benchmark, ext, README) was migrated to the canonical name.

  - **Step-size kwarg canonicalised to `h`** (api-review §3(a).1):
    - `solve_pade(; h_max)` → `solve_pade(; h)` (still required).  Added
      a docstring note that `solve_pade` is a *fixed-step* solver — `h`
      is the exact step, clamped only at the span end.
    - `lattice_dispatch_solve(; h_path)` → `lattice_dispatch_solve(; h)`;
      the internal forwarding to `path_network_solve` /
      `edge_gated_pole_field_solve` now passes a single `h`, retiring the
      apologetic `h_path`-vs-`h` docstring parenthetical.
    - `PadeTaylorAlg(; h_max)` constructor kwarg → `PadeTaylorAlg(; h)`.
  - **Iteration-cap kwargs snake-cased to `max_iter`** (api-review
    §3(a).2) — *naming only; the distinct concepts are NOT unified*:
    - `bvp_solve(; maxiter)` (both 2-arg and 3-arg overloads) →
      `bvp_solve(; max_iter)`.
    - `vector_bvp_solve(; maxiter)` → `vector_bvp_solve(; max_iter)`.
    - `dispatch_solve(; bvp_maxiter)` → `dispatch_solve(; bvp_max_iter)`;
      its internal `bvp_solve` forwarding now passes `max_iter`.
    - Untouched (already correct, distinct concepts): `max_steps`,
      `max_rescales`, `max_steps_per_target`, and
      `edge_gated_pole_field_solve`'s existing `max_iter`.

### Breaking (theoretical, pre-v1.0) — `PadeTaylorAlg` field rename

  - The `PadeTaylorAlg` struct **field** `h_max` was renamed to `h`
    (`PadeTaylorAlg{H}.h`).  Julia binds struct fields by name, so —
    unlike a kwarg — a field cannot be shimmed; this is a clean break.
    The only consumer is the in-repo `PadeTaylorCommonSolveExt`
    (`alg.h_max` → `alg.h`), updated in lockstep.  No external code is
    known to read the field directly.  The *constructor* kwarg `h_max`
    remains accepted as a deprecated alias (above), so
    `PadeTaylorAlg(; h_max = …)` still works.

### Deprecated — old kwarg names (one-release shims; beads `xds` + `0xn`)

  Each renamed kwarg accepts its OLD name as a deprecated alias that
  maps onto the new name and emits a `Base.depwarn` (the kwarg form;
  `Base.@deprecate_binding` does not cover keyword arguments).  To be
  removed in a future release.

  - `solve_pade(; h_max = …)`           → use `h`.
  - `PadeTaylorAlg(; h_max = …)`        → use `h`.
  - `lattice_dispatch_solve(; h_path = …)` → use `h`.
  - `bvp_solve(; maxiter = …)`          → use `max_iter`.
  - `vector_bvp_solve(; maxiter = …)`   → use `max_iter`.
  - `dispatch_solve(; bvp_maxiter = …)` → use `bvp_max_iter`.

  Regression guard: `test/api_kwarg_deprecation_test.jl` asserts each
  alias still produces the identical result AND emits the deprecation
  (the warn-assertion is meaningful under `Pkg.test`'s `--depwarn=yes`).

### Fixed — closed-form Painlevé constructors now exported (bead `gvz`)

  - `pii_rational`, `pii_airy`, `piv_entire` are now reachable through
    bare `using PadeTaylor` (api-review §9 row 3).  They were `export`ed
    from the `Painleve` submodule but omitted from the main module's
    selective `using .Painleve: …` import and top-level `export` block,
    so `using PadeTaylor; pii_rational(1)` threw `UndefVarError`.
    Regression-guarded by `test/painleve_closed_form_test.jl` CF.5.1
    (asserts main-module resolution, not the qualified
    `PadeTaylor.Painleve.…` path the other CF tests use).

### B5 follow-ups — all shipped

  - FFW Fig 3 (PVI phase portraits, bead `padetaylor-a1l`) and Fig 7
    (generic PVI in η/ζ/z planes, bead `padetaylor-mgx`) both closed;
    see `figures/ffw2017_fig_{3,7}.jl`, `test/ffw_fig_{3,7}_test.jl`,
    worklogs 046/047.  (Section was stale until 2026-08-23.)

## [0.1.0] — 2026-05-13

First numbered release.  All four architectural layers and five
composition tiers shipped; 1311 / 1311 tests passing.

### Added

**Tier 1 — the four-layer core (Phases Z, 1–6).**

  - `LinAlg` — SVD dispatcher with relative-accuracy guarantees on
    small singular values (LAPACK for `Float64`,
    `GenericLinearAlgebra` one-sided Jacobi for `BigFloat`,
    `Matrix{Arb} → BigFloat` shim via the `PadeTaylorArblibExt`
    extension).
  - `RobustPade.robust_pade` + `PadeApproximant` — GGT 2013
    Algorithm 2 with Chebfun's QR-reweighting (lines 111–117 of
    `padeapprox.m`).
  - `Coefficients.taylor_coefficients_1st` /
    `taylor_coefficients_2nd` — Taylor jet generation via
    `TaylorSeries.jl::Taylor1{T}` (validated for
    `T ∈ {Float64, BigFloat, Arblib.Arb}`).
  - `StepControl.step_jorba_zou` — Jorba–Zou 2005 §3.3.1 eq. 11
    truncation-error step formula (TaylorIntegration.jl-equivalent).
  - `StepControl.step_pade_root` — FW 2011 §3.1 forward-projection
    Padé-denominator-root distance heuristic.
  - `PadeStepper.pade_step!` and `pade_step_with_pade!` — one-step
    orchestrator from `(z, u, u')` to `(z + h, u(z + h), u'(z + h))`.
  - `Problems.PadeTaylorProblem` / `solve_pade` /
    `PadeTaylorSolution` / `taylor_eval` — public driver layer with
    dense output over multi-segment trajectories.

**Tier 1.5 — package extensions (Phases 7–8).**

  - `PadeTaylorCommonSolveExt` — `CommonSolve.jl` adapter exposing
    `PadeTaylorAlg <: CommonSolve.AbstractAlgorithm` with
    `init` / `step!` / `solve!` methods so PadeTaylor.jl participates
    in the wider SciML solver ecosystem.
  - `PadeTaylorArblibExt` — Arbitrary-precision SVD via
    `Arb → BigFloat → GenericLinearAlgebra.svd` (Arblib.jl has no
    native SVD; verified by source inspection).

**Tier 2 — path-network (Phase 10).**

  - `PathNetwork.path_network_solve` + `PathNetworkSolution` — FW 2011
    §3.1 five-direction wedge path-tree solver with Stage-2
    fine-grid barycentric extrapolation from stored Padé
    approximants.  Generic in `T <: AbstractFloat` for `Float64` and
    `Complex{T}`.  Optional Schwarz-reflection symmetry mode
    (`enforce_real_axis_symmetry::Bool` kwarg) for real-coefficient
    real-IC problems.

**Tier 2.5 — pole-field edge detection (Phase 12.5).**

  - `EdgeDetector.laplacian_residual` + `pole_field_mask` — 5-point
    Laplacian classifier per FW 2011 §3.2.2 (threshold on
    `log₁₀|Δu|`).

**Tier 3 — boundary-value composition (Phases 11, 12).**

  - `BVP.bvp_solve` + `BVPSolution` — Chebyshev spectral-collocation
    Newton solver for second-order analytic BVPs on a complex
    segment with Dirichlet BCs (FW 2011 §3.2; Trefethen *SMIM*
    chapters 6 & 13; Berrut–Trefethen barycentric evaluation).
    Step-norm Newton convergence (`eps(T)^(3/4)` default).
  - `Dispatcher.dispatch_solve` + `IVPSegment` / `BVPSegment` /
    `DispatcherSolution` — 1D IVP↔BVP chain composition layer per
    FW 2011 §4.4 with junction derivative-match diagnostics.
  - `LatticeDispatcher.lattice_dispatch_solve` + `LatticeSolution` —
    2D-grid composition with per-row BVP fill on smooth runs
    flanked by IVP cells (FW 2011 line 190).

**Tier 4 — exponential coordinate transforms (Phase 13).**

  - `CoordTransforms.pIII_transformed_rhs` /
    `pV_transformed_rhs` / `pIII_z_to_ζ` / `pIII_ζ_to_z` /
    `pV_z_to_ζ` / `pV_ζ_to_z` — Exponential coordinate maps
    (`z = exp(ζ/2)` for PIII, `z = exp(ζ)` for PV) that remove the
    fixed branch point at `z = 0`.  Helpers-only; integration via
    composition with `PadeTaylorProblem` + `path_network_solve`.

**Tier 5 — Riemann-sheet tracking (Phase 14).**

  - `SheetTracker.pVI_transformed_rhs` — ζ-plane RHS for the sixth
    Painlevé equation (FFW 2017 eq. 3).
  - `SheetTracker.winding_delta` / `accumulate_winding` /
    `sheet_index` — path-side Riemann-sheet bookkeeping primitives
    (signed angle change normalised to `(-π, π]`, cumulative
    winding, sheet-index assignment via `round(total / 2π)`).

**Documentation.**

  - Documenter.jl-built docs site at `docs/build/` (regenerable via
    `julia --project=docs docs/make.jl`).  Sections: Home,
    Architecture, API, Figures.  Local-only per CLAUDE.md Rule 11
    (no `deploydocs`, no CI).
  - Four Architecture Decision Records under `docs/adr/`:
    four-layer architecture (0001), BigFloat-SVD via
    `GenericLinearAlgebra` (0002), Pkg.jl weak-dep extensions (0003),
    path-network architecture (0004).
  - Figure-acceptance catalogue at `docs/figure_catalogue.md`
    covering 79 figures across FW 2011, FW 2014, FW 2015,
    RF 2014, FFW 2017.
  - 18 worklog shards under `docs/worklog/` recording frictions
    surfaced + mutation-proof procedures + algorithmic findings for
    each shipped phase.

### Headline empirical results

  - **FW 2011 Table 5.1 row 1** (long-range integration of the
    equianharmonic Weierstrass ℘-function to `z = 30`):
    **`2.13·10⁻¹⁴` rel-err** in `BigFloat`-256 — beats FW's
    published `8.34·10⁻¹⁴`.
  - **FW 2011 Fig 4.1 step (i)** (tritronquée BVP on `[-20i, +20i]`):
    `u(0)` pinned to `≤ 3.5·10⁻¹³`, `u'(0)` to `≤ 5.3·10⁻¹¹` vs
    FW eq. 4.1 reference values at `N = 240` Chebyshev nodes.
  - **FW 2011 Fig 3.1** (PI tritronquée pole field): qualitative
    PARTIAL reproduction at 25×25 over `[-4, 4]²` — 4-of-5 pole-free
    sectors recovered, conjugate symmetry verified, leading-pole
    magnitude matches Joshi–Kitaev to `≤ 10⁻³`.
  - **Phase-6 pole-bridge demo**: at `z = 1.05` (just past the
    Weierstrass-℘ pole at `z = 1`), the order-30 Padé conversion
    matches the closed-form solution to `3.45·10⁻¹⁰` while plain
    Taylor truncation of the same coefficients diverges to `2.5`
    (relative error) — a `9.86`-orders-of-magnitude gap on
    identical input.

### Cross-validation oracles

The test suite cross-validates against:

  - Mathematica's closed-form `WeierstrassP[z + c₁, {0, c₂}]`;
  - Mathematica's `NDSolve` at `WorkingPrecision = 50`;
  - `mpmath.odefun` at 40 decimal digits (Python);
  - Chebfun's `padeapprox.m` under Octave (Phase 2 oracle);
  - DMSUITE `chebdif` / `chebint` under Octave (Phase 11 BVP oracle);
  - `TaylorIntegration.jl::stepsize` + mpmath + wolframscript
    (Phase 4 three-source pin).

### Project discipline

  - 14 source modules, each under the 200-LOC ceiling (CLAUDE.md
    Rule 6).
  - All load-bearing tests mutation-proven (CLAUDE.md Rule 4):
    perturb impl, confirm RED, restore.
  - Ground-truth-before-code (CLAUDE.md Law 1) enforced via
    line-cited references in commit messages and ADRs to
    `references/markdown/<paper>/<file>.md`.
  - No GitHub CI (CLAUDE.md Rule 11) — quality gates run locally
    via `julia --project=. -e 'using Pkg; Pkg.test()'`.

### Known limitations (v1)

  - **Coordinate-transform tier (Phase 13)** ships RHS factories
    and IC round-trips for PIII / PV.  Non-uniform Stage-1 node
    placement and adaptive Padé `h` are independent follow-ons,
    not bundled.  Suitable for FFW 2017 Fig 1, 4, 5, 6 at one-step
    accuracy.
  - **Sheet-tracking tier (Phase 14)** ships the ζ-plane PVI RHS
    and post-walk winding-number primitives.  Constrained-wedge
    `PathNetwork` routing (that refuses to overstep branch cuts
    during the walk) and sheet-aware Stage-2 evaluation are
    deferred.  Suitable for FFW 2017 Fig 2, 3, 7 at one-step
    accuracy.
  - **`Polynomials.roots` for `Arb` element type** — validated only
    for `Float64` / `Complex{Float64}`.  Arb-precision Padé-root
    step deferred (friction bead `padetaylor-8pi`).
  - **GitHub-pages deployment** — explicitly out of scope per
    CLAUDE.md Rule 11.

### References

  - B. Fornberg & J. A. C. Weideman, *A numerical methodology for
    the Painlevé equations*, J. Comput. Phys. 230 (2011), 5957–5973.
  - P. Gonnet, S. Güttel & L. N. Trefethen, *Robust Padé
    Approximation via SVD*, SIAM Review 55 (2013), 101–117.
  - À. Jorba & M. Zou, *A software package for the numerical
    integration of ODEs by means of high-order Taylor methods*,
    Experimental Mathematics 14 (2005), 99–117.
  - M. Fasondini, B. Fornberg & J. A. C. Weideman, *Methods for
    the computation of the multivalued Painlevé transcendents on
    their Riemann surfaces*, J. Comput. Phys. 344 (2017), 36–50.

[Unreleased]: https://github.com/tobiasosborne/PadeTaylor.jl/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tobiasosborne/PadeTaylor.jl/releases/tag/v0.1.0
