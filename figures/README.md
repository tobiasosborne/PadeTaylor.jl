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
```

Each script prints timing + a one-line acceptance note and writes its
PNG to `figures/output/`.

## Figures

| script | source | what it shows |
|---|---|---|
| `fw2011_fig_3_1.jl` | FW 2011 Fig 3.1 (`...md:143-145`) | `\|u(z)\|` pole-field surface for PI, near-tritronquée ICs, over `[-10,10]²` |
| `fw2011_fig_3_2.jl` | FW 2011 Fig 3.2 (`...md:160-164`) | the Stage-1 path tree: `40×40` coarse grid, `h = 0.3`, rooted at the origin |

`...md` is `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md`.

## Acceptance

Per `docs/figure_catalogue.md`: Fig 3.1 is a qualitative pole-field
match (pole-free near-tritronquée region around the origin, pole
fields in the five main sectors, four of them displaced from the
origin); Fig 3.2 is a visual topology match (a connected tree rooted
at the origin reaching within `h` of every coarse-grid node, no
crossing paths). Neither figure carries a quantitative criterion in
the paper itself.
