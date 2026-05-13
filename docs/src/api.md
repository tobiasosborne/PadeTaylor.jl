# API reference

This page is generated from the in-source docstrings of each module
under `src/`.  The order mirrors `src/PadeTaylor.jl`'s `include`
sequence: the four-layer core (`LinAlg` → `RobustPade` →
`Coefficients` → `StepControl` → `PadeStepper`), then the public
driver layer (`Problems`), then the composition tiers
(`PathNetwork`, `BVP`, `Dispatcher`, `EdgeDetector`,
`LatticeDispatcher`), then the helpers for higher Painlevé equations
(`CoordTransforms`, `SheetTracker`).

For a narrative-level introduction to the four-layer architecture and
the rationale behind each module, see the [Architecture](architecture.md)
chapter.

## Umbrella module

```@autodocs
Modules = [PadeTaylor]
Order   = [:module, :type, :function]
```

## Tier 1 — the four-layer core

### `LinAlg` — SVD dispatcher

```@autodocs
Modules = [PadeTaylor.LinAlg]
Order   = [:module, :type, :function]
```

### `RobustPade` — GGT 2013 + Chebfun reweighting

```@autodocs
Modules = [PadeTaylor.RobustPade]
Order   = [:module, :type, :function]
```

### `Coefficients` — Taylor jet generation

```@autodocs
Modules = [PadeTaylor.Coefficients]
Order   = [:module, :type, :function]
```

### `StepControl` — step-size selection

```@autodocs
Modules = [PadeTaylor.StepControl]
Order   = [:module, :type, :function]
```

### `PadeStepper` — one-step orchestrator

```@autodocs
Modules = [PadeTaylor.PadeStepper]
Order   = [:module, :type, :function]
```

## Tier 1 driver — `Problems`

```@autodocs
Modules = [PadeTaylor.Problems]
Order   = [:module, :type, :function]
```

## Tier 2 — `PathNetwork`

```@autodocs
Modules = [PadeTaylor.PathNetwork]
Order   = [:module, :type, :function]
```

## Tier 3 — BVP and 1D/2D composition

### `BVP` — Chebyshev–Newton spectral collocation

```@autodocs
Modules = [PadeTaylor.BVP]
Order   = [:module, :type, :function]
```

### `Dispatcher` — 1D IVP↔BVP chain

```@autodocs
Modules = [PadeTaylor.Dispatcher]
Order   = [:module, :type, :function]
```

### `EdgeDetector` — Laplacian pole-field classifier

```@autodocs
Modules = [PadeTaylor.EdgeDetector]
Order   = [:module, :type, :function]
```

### `LatticeDispatcher` — 2D-grid composition

```@autodocs
Modules = [PadeTaylor.LatticeDispatcher]
Order   = [:module, :type, :function]
```

## Tier 4 — `CoordTransforms` (PIII, PV)

```@autodocs
Modules = [PadeTaylor.CoordTransforms]
Order   = [:module, :type, :function]
```

## Tier 5 — `SheetTracker` (PVI)

```@autodocs
Modules = [PadeTaylor.SheetTracker]
Order   = [:module, :type, :function]
```
