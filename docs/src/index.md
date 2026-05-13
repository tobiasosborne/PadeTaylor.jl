# PadeTaylor.jl

A general-purpose Taylor–Padé IVP solver for analytic ODEs `y' = f(z, y)`
with `z ∈ ℂ`, in the Fornberg–Weideman tradition.  The local solution is
expanded as a high-order Taylor series, converted to a robust Padé
rational via a singular-value-decomposition method, then evaluated on
the rescaled unit interval.  The Padé conversion provides an analytic
continuation of the local jet whose domain of validity extends *past*
the radius of convergence of the underlying Taylor series, allowing the
integrator to step across (or close to) singularities of the solution
where naïve truncation diverges.

## Status

**Research-grade, in progress.**  All four architectural layers and
five composition tiers are shipped; the package is not yet registered
in the Julia general registry.  1311 / 1311 tests passing as of the
SheetTracker PVI ζ-plane commit.

| Phase | Module / Feature | Status |
|---|---|---|
| Z | scaffold + ADRs + project discipline | ✅ shipped |
| 1 | `LinAlg.pade_svd` | ✅ shipped |
| 2 | `RobustPade.robust_pade` | ✅ shipped |
| 3 | `Coefficients.taylor_coefficients_*` | ✅ shipped |
| 4 | `StepControl.step_jorba_zou` + `step_pade_root` | ✅ shipped |
| 5 | `PadeStepper.pade_step!` | ✅ shipped |
| 6 | `Problems.PadeTaylorProblem` + `solve_pade` | ✅ shipped |
| 7 | `CommonSolveAdapter` extension (SciML) | ✅ shipped |
| 8 | `PadeTaylorArblibExt` extension (Arb arbitrary precision) | ✅ shipped |
| 9 | PI tritronquée pole-field qualitative reproduction | ✅ shipped |
| 10 | `PathNetwork.path_network_solve` (Tier 2) | ✅ shipped |
| 11 | `BVP.bvp_solve` (Chebyshev–Newton, Tier 3) | ✅ shipped |
| 12 | `Dispatcher.dispatch_solve` + `LatticeDispatcher` (Tier 3) | ✅ shipped |
| 13 | `CoordTransforms` (PIII / PV, Tier 4) | ✅ shipped |
| 14 | `SheetTracker` (PVI ζ-plane + winding, Tier 5) | ✅ shipped |

The headline empirical result for the FW Table 5.1 long-range
integration of the equianharmonic ℘-function to `z = 30` is `2.13e-14`
relative error in `BigFloat`-256 — better than the `8.34e-14` reference
value reported in Fornberg & Weideman 2011.

## Installation

PadeTaylor.jl requires Julia 1.10 or newer.  The package is not yet
registered; install from this repository:

```julia
] add https://github.com/tobiasosborne/PadeTaylor.jl
```

Or for development:

```julia
] dev https://github.com/tobiasosborne/PadeTaylor.jl
```

## Headline: Padé bridges the pole, plain Taylor diverges

The canonical Phase-6 test problem is `u'' = 6u²` with the
Fornberg–Weideman initial conditions at `z = 0`; the closed-form
solution is `u(z) = ℘(z + c₁; 0, c₂)` on the equianharmonic Weierstrass-℘
lattice with `c₁ = -1`, `c₂ = 2`.  The solution has a double pole at
`z = 1`.

```julia
using PadeTaylor

f(z, u, up) = 6u^2

# ONE local Padé built at z = 0 with h_max = 1.5 covers the lattice
# pole at z = 1 (which sits at the rescaled t = 1/1.5 ≈ 0.667, strictly
# inside the segment).
prob = PadeTaylorProblem(f, (1.071822516416917, 1.710337353176786),
                         (0.0, 1.5); order = 30)
sol  = solve_pade(prob; h_max = 1.5)

sol(0.5)    # ≈ (4.0044,    15.9643)   — before the pole
sol(0.95)   # ≈ (400.00,    15999.9)   — just before
sol(1.05)   # ≈ (400.00,   -15999.9)   — just AFTER the pole; matches ℘
sol(1.4)    # ≈ (6.2518,   -31.2317)   — further past, still tracking ℘
```

At `z = 1.05`, the Padé conversion of the order-30 Taylor jet matches
the closed-form Weierstrass-℘ to relative error `3.45·10⁻¹⁰`.  Plain
Taylor truncation of the same coefficients at the same point gives a
relative error of `2.5` (i.e. 250 %) — a `9.86`-orders-of-magnitude gap
between the two paths on identical input.

```julia
coefs = taylor_coefficients_2nd(f, 0.0, 1.071822516416917,
                                   1.710337353176786, 30)
PadeTaylor.taylor_eval(coefs, 1.05)   # ≈ 1397   (true value: 400)
```

## Going further: path-networks, BVPs, and Riemann sheets

The single-segment Padé bridge above is just Tier 0 / Phase 6.  The
higher tiers let you:

- **Tier 2 (`PathNetwork`)** — Navigate the complex `z`-plane with a
  five-direction wedge of candidate steps, building a visited-node
  tree of local Padé approximants.  A subsequent Stage-2 pass fills a
  fine output grid by barycentric evaluation of the nearest stored
  Padé.  This is the Fornberg–Weideman §3.1 long-range integration
  recipe; on FW Table 5.1 it converges to `2.13·10⁻¹⁴` rel-err.

- **Tier 3 (`BVP`, `Dispatcher`, `LatticeDispatcher`)** — Chebyshev
  spectral-collocation Newton solver for smooth-band BVPs (FW 2011
  §3.2 + Trefethen SMIM), 1D IVP↔BVP chain composition, and 2D
  lattice composition with per-row BVP fill on smooth runs flanked by
  pole fields (FW 2011 line 190 recipe).

- **Tier 4 (`CoordTransforms`)** — Exponential coordinate maps that
  remove the fixed branch point at `z = 0` for PIII and PV per
  Fasondini–Fornberg–Weideman 2017 §2.1.

- **Tier 5 (`SheetTracker`)** — PVI ζ-plane transformed RHS plus
  winding-number primitives for tracking Riemann sheet index as
  traversed paths circumnavigate the branch-point lattice at
  `ζ = 2π·i·k`.

See the [Architecture](architecture.md) chapter for the design
rationale, and the [API](api.md) reference for module-by-module
docstrings.

## Tests

```julia
] test PadeTaylor
```

or, equivalently, from the project root:

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
```

The test suite cross-validates against:

- Mathematica's closed-form `WeierstrassP[z + c₁, {0, c₂}]`;
- Mathematica's `NDSolve` at `WorkingPrecision = 50`;
- `mpmath.odefun` at 40 decimal digits (Python);
- Chebfun's `padeapprox.m` under Octave (Phase 2 oracle);
- DMSUITE `chebdif`/`chebint` under Octave (Phase 11 BVP oracle).

Per-module pinned oracles live under `external/probes/<probe-name>/`
with capture scripts that re-derive the values on demand.  Quality
gates run locally; there is no GitHub CI by design (failure-noise from
automated CI is worse than zero signal at the present development
stage).

## Repository documentation

In addition to this site:

- `RESEARCH.md` — Stage 0 deliverable: deep dive on FW 2011, GGT 2013,
  the Jorba–Zou step-size formula, the BigFloat-SVD landscape, the
  `TaylorSeries.jl::Taylor1{Arb}` empirical confirmation, and
  resolutions of all open-question points.
- `DESIGN.md` — Stage 1 deliverable: the granular execution plan.
- `CLAUDE.md` — project discipline (ground-truth-before-code, TDD,
  mutation-proof, literate programming, ≤200 LOC per module).
- `HANDOFF.md` — current state for the next agent picking up.
- `docs/adr/` — Architecture Decision Records (synthesised in
  [Architecture](architecture.md)).
- `docs/worklog/` — per-phase frictions surfaced and lessons learnt.
- `docs/figure_catalogue.md` — full 79-figure tiered acceptance
  catalogue (synthesised in [Figures](figures.md)).

The `references/` directory carries the load-bearing PDFs (FW 2011,
GGT 2013, FW 2014/15, RF 2014, FFW 2017, Jorba–Zou 2005, Mezzarobba
2019), each with a `marker_single`-converted markdown extract under
`references/markdown/<name>/` for line-cited reasoning in commits and
ADRs.

## References

- B. Fornberg & J. A. C. Weideman, *A numerical methodology for the
  Painlevé equations*, J. Comput. Phys. 230 (2011), 5957–5973.
- P. Gonnet, S. Güttel & L. N. Trefethen, *Robust Padé Approximation
  via SVD*, SIAM Review 55 (2013), 101–117.
- À. Jorba & M. Zou, *A software package for the numerical
  integration of ODEs by means of high-order Taylor methods*,
  Experimental Mathematics 14 (2005), 99–117.
- M. Fasondini, B. Fornberg & J. A. C. Weideman, *Methods for the
  computation of the multivalued Painlevé transcendents on their
  Riemann surfaces*, J. Comput. Phys. 344 (2017), 36–50.

## License

PadeTaylor.jl is licensed under the GNU Affero General Public License,
version 3 or any later version (AGPL-3.0-or-later).  See `LICENSE` for
the full text.
