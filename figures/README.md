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
julia --project=figures figures/fw2011_fig_4_4.jl
julia --project=figures figures/fw2011_fig_5_1.jl
julia --project=figures figures/fw2011_fig_5_2.jl
```

Each script prints timing + a one-line acceptance note and writes its
PNG to `figures/output/`.

`figutil.jl` is a shared `include`-d helper (lattice construction +
the `Axis3` pole-field surface renderer) used by the pole-field
surface figures; it is not a runnable script.

## Figures

| script | source | what it shows |
|---|---|---|
| `fw2011_fig_3_1.jl` | FW 2011 Fig 3.1 (`...md:143-145`) | `\|u(z)\|` pole-field surface for PI, near-tritronquée ICs, over `[-10,10]²` |
| `fw2011_fig_3_2.jl` | FW 2011 Fig 3.2 (`...md:160-164`) | the Stage-1 path tree: `40×40` coarse grid, `h = 0.3`, rooted at the origin |
| `fw2011_fig_3_3.jl` | FW 2011 Fig 3.3 (`...md:202-208`) | `log₁₀\|Δu\|` pole-field edge detector (5-point Laplacian stencil) + level-`0.001` contour, tritronquée ICs |
| `fw2011_fig_4_2.jl` | FW 2011 Fig 4.2 (`...md:231-241`) | real-axis `u(x)` curves: tronquée + near-tronquée cases with `u(0)=0`, two panels, `±√(-x/6)` leading-term branches |
| `fw2011_fig_4_3.jl` | FW 2011 Fig 4.3 (`...md:243-245`) | `\|u(z)\|` pole-field surface, NIST Handbook example `u(0)=0, u'(0)=1.8518` |
| `fw2011_fig_4_4.jl` | FW 2011 Fig 4.4 (`...md:253-255`) | `\|u(z)\|` pole-field surface, `u(0)=0, u'(0)=1.8519` — the tronquée-transition companion to Fig 4.3 |
| `fw2011_fig_5_1.jl` | FW 2011 Fig 5.1 (`...md:281-318`) | Weierstrass-℘ test problem: analytic ℘ pole lattice + the Padé-integrator path to `z=30` |
| `fw2011_fig_5_2.jl` | FW 2011 Fig 5.2 (`...md:320-326`) | `log₁₀`(rel-err) surface + accuracy/time contours over the `(order,h)` plane — the `order=30, h=0.5` justification |

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
not part of the figure scripts. Fig 5.1 is a match for the ℘ pole
lattice + integrator path; Fig 5.2 meets the catalogue's quantitative
criterion — the `(order,h)` sweep minimum lands at `(30, 0.40)`
(FW's `(30, 0.5)` within `±5/±0.1`) with min rel-err `9.1e-15`.
