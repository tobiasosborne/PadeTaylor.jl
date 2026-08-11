# `figures/` — FW-family figure reproduction

This directory reproduces the "pretty pictures" of the Painlevé
methodology papers as **runnable scripts**, one script per figure.
Each script is exposition: it states the source figure, the ODE and
initial conditions, the grid and step parameters, and writes a PNG
into `figures/output/`.

The directory is its own Julia project (`figures/Project.toml`) so
that the heavyweight plotting dependency (`CairoMakie`) never enters
the `PadeTaylor.jl` package's own `Project.toml`. `PadeTaylor` itself
is a `dev`-ed path dependency, so the scripts always exercise the
working-tree source.

## Setup (one-time)

```bash
julia --project=figures -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
```

## Running a figure

```bash
julia --project=figures figures/fw2011_fig_3_1.jl
julia --project=figures figures/fw2011_fig_3_2.jl
julia --project=figures figures/fw2011_fig_3_3.jl
julia --project=figures figures/fw2011_fig_4_2.jl
julia --project=figures figures/fw2011_fig_4_3.jl
julia --project=figures figures/fw2011_fig_4_1.jl
julia --project=figures figures/fw2011_fig_4_4.jl
julia --project=figures figures/fw2011_fig_4_7.jl
julia --project=figures figures/fw2011_fig_4_8.jl
julia --project=figures figures/fw2011_fig_5_1.jl
julia --project=figures figures/fw2011_fig_5_2.jl
julia --project=figures figures/ffw2017_fig_6.jl
julia --project=figures figures/noumi_yamada_a4_pole_field.jl
julia --project=figures figures/kkg_pi2_tritronquee_pole_field.jl
julia --project=figures figures/kkg_pi2_tritronquee_surface.jl
```

Each script prints timing + a one-line acceptance note and writes its
PNG to `figures/output/`.

`figutil.jl` is a shared `include`-d helper (lattice construction, the
`Axis3` pole-field surface renderer, and the `pole_scatter_axis`
pole-location scatter renderer) used by the figure scripts; it is not
a runnable script.

## Figures

| script | source | what it shows |
|---|---|---|
| `fw2011_fig_3_1.jl` | FW 2011 Fig 3.1 (`...md:143-145`) | `\|u(z)\|` pole-field surface for PI, near-tritronquée ICs, over `[-10,10]²` |
| `fw2011_fig_3_2.jl` | FW 2011 Fig 3.2 (`...md:160-164`) | the Stage-1 path tree: `40×40` coarse grid, `h = 0.3`, rooted at the origin |
| `fw2011_fig_3_3.jl` | FW 2011 Fig 3.3 (`...md:202-208`) | `log₁₀\|Δu\|` pole-field edge detector (5-point Laplacian stencil) + level-`0.001` contour, tritronquée ICs |
| `fw2011_fig_4_2.jl` | FW 2011 Fig 4.2 (`...md:231-241`) | real-axis `u(x)` curves: tronquée + near-tronquée cases with `u(0)=0`, two panels, `±√(-x/6)` leading-term branches |
| `fw2011_fig_4_3.jl` | FW 2011 Fig 4.3 (`...md:243-245`) | `\|u(z)\|` pole-field surface, NIST Handbook example `u(0)=0, u'(0)=1.8518` |
| `fw2011_fig_4_1.jl` | FW 2011 Fig 4.1 (`...md:214-229`) | near-tritronquée `\|u(z)\|` surface composed by FW's 3-step recipe: (i) imaginary-axis `bvp_solve`, (ii) two `edge_gated_pole_field_solve` run-outs, (iii) per-row `bvp_solve` smooth-band fill |
| `fw2011_fig_4_4.jl` | FW 2011 Fig 4.4 (`...md:253-255`) | `\|u(z)\|` pole-field surface, `u(0)=0, u'(0)=1.8519` — the tronquée-transition companion to Fig 4.3 |
| `fw2011_fig_4_7.jl` | FW 2011 Fig 4.7 (`...md:307-308`) | six PI pole-location scatters over `[-50,50]²` for six IC choices; each via `edge_gated_windowed_poles` (ADR-0034 — windowed composite for a seam-free pole set, edge-gated mask filter for a bloom-free one); panel (f) is the sector-confined tritronquée |
| `fw2011_fig_4_8.jl` | FW 2011 Fig 4.8 (`...md:312-314`) | PI pole-location scatter, `u(0)=-5, u'(0)=0`, over `[-90,30]×[-30,30]` — the sharp pattern transition near `Re(z)≈-60` |
| `fw2011_fig_5_1.jl` | FW 2011 Fig 5.1 (`...md:281-318`) | Weierstrass-℘ test problem: analytic ℘ pole lattice + the Padé-integrator path to `z=30` |
| `fw2011_fig_5_2.jl` | FW 2011 Fig 5.2 (`...md:320-326`) | `log₁₀`(rel-err) surface + accuracy/time contours over the `(order,h)` plane — the `order=30, h=0.5` justification |
| `ffw2017_fig_6.jl` | FFW 2017 Fig 6 (`references/markdown/FFW2017_painleve_riemann_surfaces_preprint/FFW2017_painleve_riemann_surfaces_preprint.md:295-297`) | generic `P_V` `(α,β,γ,δ)=(1,−1,1,−1/2)` on three Riemann sheets — nine panels (`|w(ζ)|`, `|u(z)|`, `arg u(z)` × sheets 0/1/2); the first FFW 2017 figure reproduced |
| `noumi_yamada_a4_pole_field.jl` | v0.2 PRD `n ≥ 4` acceptance (`docs/v0p2_plan.md` row V8a; pillar B §1, §5.2) | A_4⁽¹⁾ Noumi–Yamada pole field — generic α, all-nonzero IC; three panels (shared-`Q` pole scatter, Stage-1 visited tree, `log₁₀‖f‖` node field); the first vector-stack higher-Painlevé pole-field figure |
| `kkg_pi2_tritronquee_pole_field.jl` | KKG 2015 Figs 7.4/7.5 (`references/painleve_hierarchy/KapaevKleinGrava2015_PI2_tritronquee_ConstrApprox41.pdf`; v0.2 plan row V8b; pillar C §1, §4) | P_I⁽²⁾ tritronquée V₀ (t = 0) — Panel A the vector-BVP solution on the negative real axis `[-20,-2]` (2+2 (u,u') split, matches the KKG `n=2` asymptotic to `3.4e-4`), Panel B the pole-field scatter (21 poles in the `arg x ≈ 0` wedge from a path-network march; 389 visited nodes); the V8b vector-epic figure |
| `kkg_pi2_tritronquee_surface.jl` | KKG 2015 Figs 7.4/7.5 (`references/painleve_hierarchy/KapaevKleinGrava2015_PI2_tritronquee_ConstrApprox41.pdf`; v0.2 plan row V8b+; ADR-0024; pillar C §1, §4) | **the v0.2 headline figure** — the full `Re V₀` (KKG Fig 7.4) and `Im V₀` (KKG Fig 7.5) of the P_I⁽²⁾ tritronquée over the complex x-plane `[-20,20]²`, each as a 2D `:RdBu`-diverging heatmap and a 3D `Axis3` surface (2×2 panel grid). The smooth ~270° sector is the F3 triple-method majority vote (ray-fan BVP / 2D-Chebyshev + Gridap FEM Laplace); the `arg x ≈ 0` wedge is the path-network + Stage-2 fill (`\|u\|` spikes to O(10³), display-clamped to ±15 → the KKG jagged plateau); 21 wedge poles overlaid as dots; the negative-real-axis ridge is positive (`Re V₀(-15,0) = +4.48`, KKG eq. (1.3)) |

`...md` is `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md`.

## Acceptance

Per `docs/figure_catalogue.md`: Fig 3.1 is a qualitative pole-field
match (pole-free near-tritronquée region around the origin, pole
fields in the five main sectors, four of them displaced from the
origin); Fig 3.2 is a visual topology match (a connected tree rooted
at the origin reaching within `h` of every coarse-grid node, no
crossing paths); Fig 3.3 is a visual match for the edge-detector
signature (a deep flat smooth plain, sharp pole-field ridges, the
level-`0.001` contour cleanly separating the two); Fig 4.2 is a match
for the plotted real-axis curves (tronquée curves hugging the
`±√(-x/6)` branches, near-tronquée curves carrying real-axis poles);
Fig 4.3 / 4.4 are visual pole-field matches whose near-identical pair
brackets the tronquée transition. The fine quantitative pole-count
pin for Fig 4.3/4.4 is a test-suite obligation (bead `padetaylor-p3l`),
not part of the figure scripts. Fig 4.1 is a visual `|u(z)|`-surface
match for FW's 3-step BVP+IVP+BVP composition, cross-validated by
`test/fw_fig_41_test.jl` `FF.2.*` (the step-(ii) IVP run-out and the
step-(i) BVP spine agree to ≤1e-6 — independent methods). Fig 4.7 is a six-panel visual
pole-location match; both it and Fig 4.8 are driven by
`edge_gated_windowed_poles` (ADR-0034), which composes the bounded-window
composite (seam-free: pole-set seed-invariance 77%→99.4%, gated in-suite by
`test/field_seam_test.jl`) with the edge-gated pole-field mask (bloom-free:
off-wedge poles 4475→0 on panel (f), the sector-confined tritronquée).
Fig 4.8 is a visual match for the pattern transition near `Re(z)≈-60`.
Fig 5.1 is a match for the ℘ pole
lattice + integrator path; Fig 5.2 meets the catalogue's quantitative
criterion — the `(order,h)` sweep minimum lands at `(30, 0.40)`
(FW's `(30, 0.5)` within `±5/±0.1`) with min rel-err `9.1e-15`.

`noumi_yamada_a4_pole_field.jl` and `kkg_pi2_tritronquee_pole_field.jl`
are the v0.2 vector-epic figures; there is no published
machine-readable reference dataset for either, so their quantitative
obligations live in the matching `test/*_figure_test.jl` files
(`NYF.*`, `PI2F.*`) which reproduce the figure's exact computation via
the shared Makie-free helper kernel and assert genuine invariants.
For `kkg_pi2_tritronquee_pole_field.jl`: Panel A (the vector BVP) is
the proven quantitative deliverable (probe `padetaylor-0ln.29`
`RECIPE.md`; BVP residual `1.5e-11`, KKG `n=2` asymptotic match
`3.4e-4`); Panel B (the pole-field march) succeeds as shipped
(`h = 0.1`, 21 poles in the `arg x ≈ 0` wedge) — the wedge step `h`
is a hand-tuned v1 corner (the principled shared-`Q`-root wedge
criterion is deferred bead `padetaylor-0ln.23`).

`kkg_pi2_tritronquee_surface.jl` is the v0.2 headline figure: the full
KKG Fig 7.4/7.5 `Re V₀`/`Im V₀` surface over `[-20,20]²`, computed by
the F3 triple-method kernel `_kkg_pi2_surface_helpers.jl` (bead
`padetaylor-0ln.33`) — itself validated by `test_kkg_pi2_surface.jl`
(17888 assertions GREEN). F4 is the render only; its acceptance is the
visual KKG-7.4/7.5 match — a smooth ~270° sector with a positive
negative-real-axis ridge, the jagged `arg x ≈ 0` pole wedge, and `Im V₀`
odd under `y ↦ -y` (Schwarz symmetry). The four F3 v1 corners
(inner-arc asymptotic seed, wedge Stage-2 extrapolation, hand-tuned
wedge `h = 0.1`, ±3° Stokes-strip masking) carry through to this figure
and are documented in its header docstring.
