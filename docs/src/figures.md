# Figure-acceptance catalogue

PadeTaylor.jl is verified against the literature primarily through
**figure reproduction** — the canonical Painlevé methodology papers
(Fornberg & Weideman 2011, FW 2014, FW 2015, Reeger & Fornberg 2014,
Fasondini–Fornberg–Weideman 2017) carry roughly 80 figures showing
pole fields, tronquée curves, pole-counting diagrams, and Riemann
surfaces.  Each is potentially a regression target.

This chapter summarises the tiered figure-acceptance plan.  The
**full per-figure catalogue** with line-cited references for all 79
figures lives at `docs/figure_catalogue.md` in the repository — when
a new tier ships, the corresponding rows graduate from "target" to
"shipped" and a worklog shard records the pinning procedure used.

## Tier definitions

| tier | algorithmic prerequisite                                                          | what it unlocks                                                                                                          |
| ---- | --------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------ |
| T0   | Single-segment Padé (Phase 6)                                                     | The Padé-vs-Taylor pole-bridge demo at one `z` step.                                                                     |
| T1   | Multi-segment fixed-`h` Padé (no path-network)                                    | Test ODEs without complex-plane navigation.                                                                              |
| T2   | **Path-network** (Phase 10) + edge detector (Phase 12.5)                          | 2D pole-field plots, IC-plane survey diagrams, pole-counting diagrams, FW Table 5.1 long-range.                          |
| T3   | T2 + **BVP solver** (Phase 11) + dispatchers (Phase 12 v1 + v2)                   | Tronquée solutions with adjacent pole-free sectors; near-tritronquée with smooth imaginary band; hybrid IVP + BVP fills. |
| T4   | T3 + **exponential coordinate transform** for PIII / PV                           | Multi-sheeted transcendents with branch points at fixed locations; single-sheet Riemann-surface views.                   |
| T5   | T4 + **sheet tracking + branch-cut routing**                                      | PVI multi-sheet views; circumambulation around branch points; phase-portrait Riemann surfaces.                           |

## Acceptance discipline

Default acceptance per figure:

  - **Pole-location plots:** all poles within the domain match FW's
    pinned locations to ≤`1e-6` absolute (`Float64`) or ≤`1e-13`
    (`BigFloat`-256).  Count, residue sign, and topology preserved.
  - **Real-axis solution profiles:** rel-err vs the three-source
    oracle (closed-form ≡ Mathematica `NDSolve` ≡ `mpmath` where
    available) ≤`1e-9` (`Float64`) or ≤`1e-13` (`BigFloat`-256).
  - **Pole-counting diagrams:** count agrees with FW's diagram cells
    over a coarse IC sample (≥20×20 grid points across the published
    window).
  - **Error / parameter surfaces** (FW 2011 Fig 5.2 type): qualitative
    shape — favourable region, monotonicity along axes — matches FW;
    rel-err of the surface minimum within one order of magnitude of
    FW's reported value.

Exact criteria per figure are stated in the per-figure tables in
`docs/figure_catalogue.md`.

## Tier-by-tier status

| tier | infrastructure status | example pinned figures | next graduation work |
|------|------------------------|------------------------|----------------------|
| **T0** | ✅ shipped (Phase 6) | Phase-6 ℘-function pole-bridge demo | — |
| **T1** | ✅ shipped (Phase 6 stepper) | FW 2011 Fig 2.1, Fig 5.2 (target) | per-figure pinning |
| **T2** | ✅ shipped (Phase 10 PathNetwork + Phase 12.5 EdgeDetector) | FW 2011 Fig 5.1 (FW Table 5.1 ≤`2.13e-14` `BF`-256, beats FW's `8.34e-14`); Fig 3.1 + Fig 3.2 rendered in `figures/` (worklog 022) | quantitative per-figure pinning across FW 2011 §3, FW 2014 + 2015, RF 2014 |
| **T3** | ✅ shipped (Phase 11 BVP + Phase 12 v1 1D Dispatcher + Phase 12 v2 LatticeDispatcher) | FW 2011 Fig 4.1 step-(i) BVP pin (`u(0) ≤ 3.5e-13`, `u'(0) ≤ 5.3e-11`, see `test/fw_fig_41_test.jl`) | Fig 4.1 steps (ii) + (iii); FW 2014 tronquée figures; FFW 2017 hybrid IVP + BVP figures |
| **T4** | 🟡 PARTIAL — `CoordTransforms` shipped (Phase 13, `src/CoordTransforms.jl` exposes PIII / PV RHS factories + IC round-trips, end-to-end direct-vs-transformed agreement ≤`1e-10` per worklog 017) | — | non-uniform Stage-1 nodes; adaptive Padé `h`; per-figure pinning of FFW 2017 Fig 1, 4, 5, 6 |
| **T5** | 🟡 PARTIAL — `SheetTracker` shipped (Phase 14, `src/SheetTracker.jl` exposes PVI ζ-plane RHS + winding primitives, end-to-end direct-vs-transformed agreement ≤`1e-10` per worklog 018) | — | η-plane PVI eq.; constrained-wedge `PathNetwork` routing that enforces branch-cut avoidance during the walk; per-figure pinning of FFW 2017 Fig 2, 3, 7 |

## Headline pinned figures

The figures pinned with concrete quantitative acceptance in the test
suite to date:

### FW 2011 Fig 5.1 — `℘`-function path to `z = 30`

The headline long-range test.  Path-network walk from `z = 0` to
`z = 30` on the equianharmonic Weierstrass lattice; FW report
`8.34e-14` relative absolute error at `z = 30` with `BF`-256 and
order-30 Padé.  PadeTaylor.jl converges to **`2.13e-14`** on the same
problem (worklog 008).  Lives in `test/pathnetwork_test.jl` under
PN.2.2.  `figures/fw2011_fig_5_1.jl` renders it: the analytic ℘ pole
lattice (large dots) and the integrator path threading the low-`|u|`
channel between the pole rows (worklog 025).

### FW 2011 Fig 5.2 — the `(order, h)` accuracy/cost trade-off

`figures/fw2011_fig_5_2.jl` sweeps `order ∈ 4:2:50` against
`h ∈ 0.05:0.05:1.0` — 480 ℘ walks to `z = 30` — and draws the
`log₁₀`(rel-err) surface plus smoothed accuracy/compute-time contours.
The sweep minimum lands at `(order, h) = (30, 0.40)` with min rel-err
`9.1e-15`, reproducing FW's central claim that `order = 30, h = 0.5`
is the favourable choice (worklog 025).

### FW 2011 Fig 4.1 step (i) — tritronquée BVP pin

Chebyshev–Newton BVP for `u'' = 6u² + z` on `[-20i, +20i]` with
leading-term `u(z) = -√(-z/6)` Dirichlet BCs.  At `N = 240` the
solver pins `u(0) ≤ 3.5e-13` and `u'(0) ≤ 5.3e-11` vs the FW eq. 4.1
reference values — well under the 1e-10 spec (worklog 016).  Lives
in `test/fw_fig_41_test.jl`.

### FW 2011 Fig 3.1, Fig 3.2, Fig 3.3 — PI pole field, path tree, edge detector

Reproduced as runnable scripts under `figures/` (worklogs 022, 023).
`figures/fw2011_fig_3_1.jl` renders the `|u(z)|` pole-field surface on
a 121×121 lattice over `[-10,10]²`; the calm-near-origin signature is
confirmed numerically by a 3.1× median-`|u|` ratio between the central
disc `|z|≤2` and the annulus `6≤|z|≤10`.  `figures/fw2011_fig_3_2.jl`
renders the FW §3.1 Stage-1 path tree on the exact 40×40 coarse grid
at `h = 0.3` — a connected, non-crossing tree rooted at the origin,
drawn from the `PathNetworkSolution.visited_parent` edge set added in
worklog 022.  `figures/fw2011_fig_3_3.jl` renders the FW §3.2.2
pole-field edge detector: the `log₁₀|Δu|` surface from the 5-point
Laplacian stencil (`EdgeDetector.laplacian_residual`) with the
level-`0.001` contour overlaid, for the precise tritronquée ICs
(eq. 4.1) over `x∈[-4,8], y∈[-6,6]` — a deep flat smooth plain, sharp
pole-field ridges, the contour separating the two.

### FW 2011 Fig 4.2, Fig 4.3, Fig 4.4 — the `u(0) = 0` tronquée family

Reproduced as runnable scripts under `figures/` (worklog 024).
`figures/fw2011_fig_4_2.jl` renders the two-panel real-axis `u(x)`
curves — the tronquée cases hugging the `±√(-x/6)` leading-term
branches, the near-tronquée `1.8518` / `1.8519` curves carrying
real-axis poles — computed (as FW does, md:233) via
`path_network_solve` on a 1-D real-axis target grid.
`figures/fw2011_fig_4_3.jl` and `figures/fw2011_fig_4_4.jl` render the
NIST-Handbook `|u(z)|` pole-field surfaces over `[-10,10]²` for
`u'(0) = 1.8518` and `1.8519`, with FW's thick real-axis line and
origin cross-line overlaid; the two are near-identical, bracketing
the tronquée transition. Their quantitative pole-count pin is a
follow-up test (bead `padetaylor-p3l`).

The quantitative regression-test companion to Fig 3.1 is the Phase 9
pin: the pole field at 25×25 over `[-4, 4]²` with 4-of-5 pole-free
sectors recovered, conjugate symmetry verified, and leading-pole
magnitude matching Joshi–Kitaev to ≤`1e-3` (worklog 012).  Lives in
`test/phase9_tritronquee_test.jl`.

## Provenance notes

A few subtleties carried in `docs/figure_catalogue.md §7` that are
worth foregrounding:

  - **Wedge half-angle convention.**  FW 2011 lines 158–159 specifies
    `±22.5°`, `±45°`; RF 2014 lines 148–153 specifies `±15°`, `±30°`.
    PadeTaylor.jl exposes this as the `wedge_angles` keyword of
    `path_network_solve` and defaults to FW 2011's
    `[-π/4, -π/8, 0, π/8, π/4]`.

  - **Edge-detector level `0.001`.**  Calibrated empirically for PI
    on an `h = 0.5` grid (FW 2011 line 208).  Other equations or
    coarser/finer grids may need re-calibration; the threshold is a
    user-tunable parameter of `pole_field_mask`.

  - **Visual-match acceptance** for figures without a quantitative
    oracle (most FW 2014 figures) is interpreted as reproducing FW's
    rendered figure to within visual indistinguishability at print
    scale.  A more rigorous criterion (e.g., pixel-space RMSE after
    re-rendering) is on the v2 backlog.

  - **`BigFloat`-256 vs `Float64`.**  Required only for FW Table 5.1
    long-range (Fig 5.1) and FFW 2017 condition-number-blow-up
    figures (Fig 5).  Tier 2/3 figures otherwise pass at `Float64`.

## Pointers

  - **Full per-figure catalogue:** `docs/figure_catalogue.md` (in
    the repository) — 79 figures across FW 2011, FW 2014, FW 2015,
    RF 2014, FFW 2017, each with line-cited references to the source
    paper's markdown extract.
  - **Per-test pinning files:** `test/pathnetwork_test.jl` (FW Table
    5.1), `test/fw_fig_41_test.jl` (Fig 4.1), `test/phase9_tritronquee_test.jl`
    (Fig 3.1 PARTIAL).
  - **Reference markdown extracts:** `references/markdown/<paper>/`
    (`marker_single`-converted from PDF), used for line-cited
    reasoning in worklogs, commit messages, and ADRs.
