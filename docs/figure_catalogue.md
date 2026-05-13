# Figure-acceptance catalogue — FW-family "pretty pictures" as regression targets

> Tracks bead `padetaylor-rgp`. Built from the four literature
> subagent reports in this session (FW2011, FW2014, FW2015,
> RF2014+FFW2017); cross-checked against the existing markdown
> extractions in `references/markdown/<paper>/`. This is a living
> document — when a new tier ships, the corresponding rows graduate
> from "target" to "shipped" and a worklog shard records the pinning
> procedure used.

## Tier definitions

| tier | algorithmic prerequisite | what it unlocks |
|---|---|---|
| **T0** | Single-segment Padé (shipped Phase 6) | The Padé-vs-Taylor pole-bridge demo at one z step. |
| **T1** | Multi-segment fixed-h Padé (no path-network) | Test ODEs without complex-plane navigation. |
| **T2** | **Path-network** (`padetaylor-1jf`) — 5-direction Stage 1 tree + Stage 2 fine grid + (optional) edge detector `padetaylor-c2p` | 2D pole-field plots; IC-plane survey diagrams; pole-counting diagrams; FW Table 5.1 long-range. |
| **T3** | T2 + **BVP solver** (`padetaylor-???` — the new bead just filed) for smooth-region segments | Tronquée solutions with adjacent pole-free sectors; near-tritronquée with smooth imaginary band; FFW-2017-style hybrid IVP+BVP. |
| **T4** | T3 + **exponential coordinate transform** (z→ζ for PIII/PV) | Multi-sheeted transcendents with branch points at fixed locations; single-sheet Riemann-surface views. |
| **T5** | T4 + **sheet tracking + branch-cut routing** | PVI multi-sheet views; circumambulation around branch points; phase-portrait Riemann surfaces. |

## Acceptance discipline

Default acceptance per figure:
  - **Pole-location plots:** all poles within domain match FW's pinned locations to ≤1e-6 absolute (Float64) or ≤1e-13 (BF-256). Count + residue sign + topology preserved.
  - **Real-axis solution profiles:** rel-err vs three-source oracle (closed-form ≡ NDSolve ≡ mpmath where available) ≤1e-9 (Float64) or ≤1e-13 (BF-256).
  - **Pole-counting diagrams:** count agrees with FW's diagram cells over a coarse IC sample (≥20×20 grid points across the published window).
  - **Error / parameter surfaces (FW 2011 Fig 5.2 type):** qualitative shape — favourable region, monotonicity along axes — matches FW; rel-err of the surface minimum within 1 order of magnitude of FW's reported value.

Exact criteria per figure are stated in the rightmost column. "Defer" = no acceptance criterion yet; record what we want from FW's caption and revisit when its tier ships.

## §1. FW 2011 — `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md`

| fig | lines | tier | what it shows | acceptance |
|---|---|---|---|---|
| 1.1 | 42–43 | — | Schematic cartoon of PI pole-field sectors | Excluded (non-numerical). |
| 2.1 | 131–135 | T1 | Taylor vs Padé on `y' = t² + y²`, near first pole | Rel-err at z=z_pole−ε ≤1e-9 for Padé arm; Taylor arm rel-err > 1e-1. (Already shipped pattern in Phase 6 test 6.1.5.) |
| 3.1 | 143–145 | T2 | PI pole field on 20×20 square, near-tritronquée ICs | Pole count + location vs RF2014-style verification; visual match. **PARTIAL — Phase 9 ships qualitative reproduction at 25×25 [-4,4]² (4-of-5 pole-free sectors, conjugate symmetry, leading-pole magnitude); FW's full 161² Fig 3.1 acceptance deferred to Phase 12 v2 (`padetaylor-k31`).** |
| 3.2 | 160–164 | T2 | PI integration path tree, h=0.3, 40×40 coarse grid | Path topology matches FW's tree visually (no quantitative criterion in paper). |
| 3.3 | 197–208 | T2+edge | log₁₀\|Δu\| for pole-field edge detection; level 0.001 | Edge-detector bitmap on 20×20 lattice matches FW's level-0.001 contour to ≤2 lattice cells off. |
| 4.1 | 218–222 | **T3** | Near-tritronquée: two pole fields + smooth imaginary-axis band; BVP solver path shown | u(0), u'(0), u(20i), u'(20i) recovered to ≤1e-10 abs; complete pole field's pole locations agree with FW Table 5.1 reference values. |
| 4.2 | 239–242 | T2 | Real-axis tronquée + near-tronquée curves, u(0)=0 | u(z) on real axis matches FW's plotted curve at sample points to ≤1e-9 rel-err. |
| 4.3 | 243–245 | T2 | PI complex pole field, u(0)=0, u'(0)=1.8518 | Pole count in 20×20 domain matches FW; pole positions ≤1e-6 abs. |
| 4.4 | 253–255 | T2 | PI pole field, u'(0)=1.8519 (adjacent tronquée) | Pole count differs by exactly 1 from Fig 4.3 (sharp transition between tronquée values). |
| 4.5 | 257–261 | **T3** | IC plane (u(0), u'(0)) with all tronquée ICs; dot = tritronquée | Tronquée curve sampled at ≥50 IC points; tritronquée dot matches FW Table 5.1's (u₀, u'₀) to ≤1e-10. |
| 4.6 | 267–269 | **T3** | Smooth-band BVP demo over z ∈ [−18, −14] | BVP solution agrees with `√(−z/6)` asymptotic + Newton correction to ≤1e-8 abs at endpoints; matches IVP-supplied derivatives to ≤1e-7 (FW's stated tol). |
| 4.7 a–f | 307–308 | T2 | Six 50×50 PI pole fields for 6 IC choices | Per-panel pole counts match FW; pole locations agree at ≤1e-6 abs. |
| 4.8 | 312–314 | T2 | PI pole field with sharp pattern transition at Re(z) ≈ −60 | Transition Re-coordinate matches FW to ≤0.2 (visual). |
| 5.1 | 316–318 | T2 | ℘-function path to z=30 with stated abs-err 8.34e-14 | **FW Table 5.1 reproduction** — the `padetaylor-8cr` headline target. Rel-err at z=30 ≤1e-13 (BF-256) per FW. |
| 5.2 | 320–326 | T2 | Error surface vs (order, h); favourable region around (30, 0.5) | Surface minimum location matches FW's (30, 0.5) ±5 / ±0.1; min rel-err within 1 order of magnitude of FW. |

## §2. FW 2014 — `references/markdown/FW2014_second_PII_exploration_FoCM14/FW2014_second_PII_exploration_FoCM14.md`

Per the FW2014 subagent: no algorithmic change vs FW 2011 in the stepper. **No quantitative acceptance table** comparable to FW 2011 Table 5.1 — validation is via figure-pattern reproduction only. Tier-2 for most pole fields; Tier-3 for tronquée-with-pole-free-sectors figures.

| fig | lines | tier | what it shows | acceptance |
|---|---|---|---|---|
| 1 | 123–127 | — | Rational solutions u₁–u₆ poles+zeros | Excluded (closed-form; not via solver). |
| 2 | 149–152 | — | Airy solutions θ=0 | Excluded (closed-form). |
| 3 | 173–176 | T2 | Airy α=5/2 pole dynamics across θ | Pole orbits match FW visually at 6 θ samples. |
| 4 | 187–190 | T2 | Pole dynamics α ∈ [5/2, 7/2] | Pole-count transitions at α=3 match FW. |
| 5 | 199–202 | T2 | Pole-counting diagrams in (u(0), u'(0))-plane | Count cells over IC plane match FW; n⁺ / n⁻ curves traced to ≤1 cell width. |
| 6 | 224–228 | T2 | Extended Fig 5 with rational at infinity | As Fig 5 + dash-dot lines at finite-α projection. |
| 7 | 244–246 | T2 | Airy curves over half-integer α pole-counting | Airy parabola u'(0)=u(0)² recovered to ≤1e-6. |
| 8 | 248–250 | T2 | Airy curves at α ∈ {5/2, 7/2, 9/2, 11/2} | Curves match FW. |
| 9 | 260–262 | T2 | HM points on pole-counting (α ∈ {0, 0.3, 0.495, 0.5}) | HM ICs match Hastings-McLeod literature values to ≤1e-8. |
| 10 | 286–288 | **T3** | pHM + sHM solutions: real-axis + pole fields | Pole-free real-axis verified by IVP; off-axis pole fields traced via path-network. |
| 11 | 300–302 | **T3** | Ablowitz-Segur pole-counting + k labels | AS curve at α=0 + k=±1 critical values match FW. |
| 12 | 334–336 | **T3** | AS near k=1 critical | Pole-pattern change at k=1.001 vs k=0.999 visible. |
| 13 | 416–418 | T2 | PII α=1/3 real-axis for 4 k values | Profiles match FW. |
| 14 | 436–438 | T2 | qAS pole-counting α ∈ {2/3, 1, 2, 3} | n⁻ / n⁺ count cells match FW. |
| 15 | 448–450 | T2 | qHM real-axis pole positions vs α | Pole positions on real axis match FW to ≤1e-6. |
| 16 | 458–460 | **T3** | Extended Fig 10 with qHM | As Fig 10 + qHM cases. |
| 17 | 470–472 | T2 | sHM/qHM pole dynamics α ∈ [5/2, 7/2] | Pole orbits match FW. |
| 18 | 484–486 | T2 | Nine PII real-axis solutions, α=2/3 | Profiles match FW. |
| 19 | 498–500 | T2 | k=0 solutions α ∈ [0, 2] | Real-axis profiles match FW. |
| 20 | 514–516 | T2 | k=0 pole fields α ∈ [2, 3] | Pole patterns match FW. |
| 21 | 528–530 | T2 | Pole-counting α=3 edge↔interior transitions | Transition lines match FW. |
| 22 | 547–549 | T2 | Pole-counting α=3/2 | As Fig 21. |
| 23 | 551–553 | T2 | Bäcklund-related solution pair | Profiles match FW. |
| 24 | 561–563 | T2 | Pole-counting α=0 / α=1/2 with mapped points | Mapping eqs.(24)–(25) verified at sample ICs. |
| 25 | 573–575 | **T3** | Six PII α=3/2 pole fields near HM/Airy gap | As Fig 10. |

## §3. FW 2015 — `references/markdown/FW2015_imaginary_PII_PhysicaD309/FW2015_imaginary_PII_PhysicaD309.md`

Per the FW2015 subagent: ImPII real axis is pole-free (residues ±i incompatible with real u on R), so real-axis paths are unobstructed; tronquée families have asymmetric pole sectors. **No new acceptance table.** Initialization is via asymptotic series at |x|≈7–10 — implementation note for path-network ICs.

| fig | lines | tier | what it shows | acceptance |
|---|---|---|---|---|
| Graphical abstract | 29 | T2 | Cover-art pole field | Defer. |
| 1 | 106–108 | — | Schematic asymptotic regimes | Excluded. |
| 2 | 110–112 | T2 | Phase portrait of convergence class vs (u(0),u'(0)) | Class assignment over IC sample matches FW. |
| 3 | 136–138 | T2 | Tronquée b pole field + real-axis profile, α=0 | Pole positions match FW; real-axis profile matches `+√(x/2)` asymptote at large x. |
| 4 | 140–142 | T2 | Six marked solutions for α=1/2, with pole fields | Per-panel pole counts + locations match FW. |
| 5 | 144–146 | T2 | Pole switching across interpolation between (a) and (b) | Transition at t-midpoint visible. |
| 6 | 181–183 | T2 | d− contour over (u(0), u'(0)) plane | Contour map matches FW visually. |
| 7 | 185–187 | **T3** | Tronquée family (a) u→0 as x→−∞ | Pole-free sector recovered via BVP; real-axis profile matches asymptote. |
| 8 | 224–226 | **T3** | Tronquée family (b) u~+√(x/2) as x→+∞ | As Fig 7. |
| 9 | 228–230 | **T3** | Tronquée family (c) u~−√(x/2) as x→+∞ | As Fig 7. |
| 10 | 251–253 | T2 | AS curve overlaid on d− contour, k=0 minimum identified | k=0 location identified to ≤1e-6 by 1D min of d−. |
| 11 | 279–281 | **T3** | AS k=0 solutions: widest pole-free real-axis neighborhood | Pole-free band radius matches FW. |
| 12 | 283–285 | T2 | Phase portraits with k markers | As Fig 2 + k labels. |
| 13 | 287–289 | T2 | AS real-axis profiles for various k | Profiles match FW at k samples. |
| 14 | 291–293 | T2 | Pole diagrams of AS solutions | Pole positions match FW. |
| 15 | 307–309 | T2 | Pole dynamics as u(0) varies between AS solutions | Pole-orbit topology matches FW. |
| A.16 | 321–323 | — | Im s(z) Borel-Ritt summation | Excluded (not via solver). |

## §4. RF 2014 — `references/markdown/ReegerFornberg2014_PIV_fundamental_domain_PhysicaD280/ReegerFornberg2014_PIV_fundamental_domain_PhysicaD280.md`

Per the RF2014 subagent: codifies FW 2011's 5-direction algorithm explicitly (at lines 148–153 with ±15°, ±30° wedge — **reconcile with FW 2011's ±22.5°, ±45°** in the unified spec). Introduces pole-counting diagrams systematically.

| fig | lines | tier | what it shows | acceptance |
|---|---|---|---|---|
| 1 | 119–121 | — | Weyl chambers in (α, √−2β)-plane | Excluded (parameter-space schematic; symbolic). |
| 2 | 155–157 | **T3** | Closed-form PIV: poles + Float64 err + 25-digit err | Float64-err map matches RF's; BF-256 err map shows ≥10 orders gain over Float64. |
| 3 | 190–192 | T2 | PIV pole-counting α=0, β=0 | Count cells match RF over (u(0), u'(0)) ∈ [-2, 2]². |
| 4 | 194–196 | — | Legend | Excluded. |
| 5 | 206–208 | — | 8-sectors schematic | Excluded. |
| 6 | 220–222 | T2 | Pole-counting at 3 boundary points | As Fig 3 × 3. |
| 7 | 254–256 | **T3** | PIV with adjacent pole-free sectors | Pole-free sector geometry matches RF; tronquée IC u₀ identified to ≤1e-8. |
| 8 | 258–260 | T2 | Pole-counting interior fundamental domain | As Fig 3. |
| 9 | 274–276 | T2 | Pole-counting positive R-axis α=0, β>0 | Count cells match RF. |
| 10 | 278–280 | T2 | Pole-counting at 6 exterior points | As Fig 3 × 6. |
| 11 | 290–292 | T2 | Pole count in (α, β)-parameter plane | Count surfaces match RF. |
| 12 | 294–296 | T2 | Zeros + poles for various integer α | Okamoto-polynomial alignment verified. |
| 13 | 302–304 | T2 | α=6, β=0 + β=±0.1 perturbations | Pole structure perturbations match RF. |
| 14 | 318–320 | T2 | Solutions normal to parabola β=−2(α−2)² | Pole-pattern change at crossing match RF. |
| 15 | 350–352 | **T3** | PIV adjacent pole-free sectors for α=0, β=−2 | As Fig 7 + dual IC. |

## §5. FFW 2017 — `references/markdown/FFW2017_painleve_riemann_surfaces_preprint/FFW2017_painleve_riemann_surfaces_preprint.md`

Per the RF2014+FFW2017 subagent: introduces five new architectural moves — exp coords, non-uniform Stage 1 nodes, adaptive Padé h, circumambulation, IVP+BVP hybrid, sheet tracking. **All seven figures are Tier 4 or Tier 5.**

| fig | lines | tier | what it shows | acceptance |
|---|---|---|---|---|
| 1 | 98–104 | **T4** | PIII non-uniform Stage 1 + adaptive Padé + 3 sheets | Pole-pattern on each of 3 sheets matches FFW; non-uniform node density correlates with pole density. |
| 2 | 193–195 | **T5** | PVI in η-plane + circumambulating in ζ-plane | Sheet-2 result matches FFW's circumambulation output; branch-cut crossing logged correctly. |
| 3 | 197–199 | **T5** | PVI phase portrait revealing tronquée | Pole-free sector matches FFW; phase-portrait winding-number consistent across sheet boundaries. |
| 4 | 234–236 | **T4** | PV tronquée + phase portrait, 3 sheets | As Fig 1 for PV. |
| 5 | 238–240 | **T4** | PIII tronquée + BVP boundary points + cond-number | PFS-exterior + BVP-interior hybrid matches FFW; cond-number heat map shows exponential growth in pole-free sector. |
| 6 | 295–297 | **T4** | PV generic 3 sheets spiraling | Pole spiral matches FFW on each sheet. |
| 7 | 299–301 | **T5** | PVI generic 3 sheets | As Fig 2 for generic case. |

## §6. Tier-aggregated work plan

| tier | figures covered | required modules | open beads |
|---|---|---|---|
| T0 | (Phase 6 ℘-function bridge) | shipped | — |
| T1 | FW2011 Fig 2.1, 5.2 | shipped | — |
| **T2** | FW2011 Fig 3.1–3.3, 4.2–4.5, 4.7–4.8, 5.1; most FW2014 + FW2015 + RF2014 | `padetaylor-1jf` PathNetwork + `padetaylor-c2p` EdgeDetector (optional) | `padetaylor-1jf`, `padetaylor-c2p` |
| **T3** | FW2011 Fig 4.1, 4.5, 4.6; FW2014 Fig 10, 11, 12, 16, 25; FW2015 Fig 7, 8, 9, 11; RF2014 Fig 2, 7, 15 | T2 + BVP module + Dispatcher | (BVP bead just filed); a Dispatcher bead to be added |
| **T4** | FFW2017 Fig 1, 4, 5, 6 | T3 + ExpCoords transform + non-uniform node placement + adaptive Padé h | New beads, P2. |
| **T5** | FFW2017 Fig 2, 3, 7 | T4 + SheetTracker + branch-cut routing | New beads, P2. |

## §7. Provenance notes

- ±22.5° vs ±15°/30° wedge: FW 2011 lines 158–159 says **±22.5°, ±45°**; RF 2014 lines 148–153 says **±15°, ±30°**. The unified spec (next deliverable) needs to pick one or expose as parameter. Default: FW 2011's 22.5°/45° per the original methodology.
- Edge-detector level 0.001: empirical for PI on h=0.5 grid (FW 2011 line 208). Needs recalibration for other equations / grids — flag as a user-tunable parameter.
- "Visual match" acceptance: when no quantitative oracle exists (most FW 2014 figures), accept by reproducing FW's `_page_N_Figure_M.jpeg` to within visual indistinguishability at print scale. A more rigorous criterion (e.g., RMSE in pixel space after re-rendering) is deferred to a v2 visual-regression bead.
- BF-256 vs Float64: required only for FW Table 5.1 long-range (Fig 5.1) and FFW 2017 condition-number-blow-up figures (Fig 5). Tier 2/3 figures otherwise pass at Float64.

## §8. Pointers

- Existing tests at `test/problems_test.jl` already demonstrate Tier 0+1 acceptance (the Phase 6 pole-bridge demo).
- Live beads: `padetaylor-8cr` (umbrella P0; this catalogue is its concrete scope), `padetaylor-1jf` (path-network), `padetaylor-c2p` (edge detector), `padetaylor-rgp` (this catalogue), plus the BVP-module bead.
- Ground truth: every figure references its `references/markdown/<paper>/<file>.md:<lines>` so any future agent can re-read the source.
