"""
    PadeTaylorArblibExt

Package-extension wiring `Arblib.Arb` / `Arblib.Acb` element types into
`PadeTaylor.LinAlg.pade_svd` per ADR-0003 + ADR-0002.  Loaded
automatically when both `PadeTaylor` and `Arblib` are present:

```julia
using PadeTaylor, Arblib

A = Arblib.Arb[1 2; 3 4]
U, S, Vt = PadeTaylor.LinAlg.pade_svd(A; full = true)
```

## Design (ADR-0003 §"PadeTaylorArblibExt.jl")

`Arblib.jl` does NOT ship an SVD primitive (verified by source
inspection — `RESEARCH.md §5.1`).  We route through `BigFloat`:

  1. Convert `Matrix{Arb}` → `Matrix{BigFloat}` via the
     `BigFloat(::Arb)` conversion that `Arblib.jl` already provides.
     **Mid-point rounding; radius discarded** (ADR-0002 caveat).
  2. Dispatch to `GenericLinearAlgebra.svd` (one-sided Jacobi
     Demmel-Veselić) for relative-accuracy guarantees on small SVs.
  3. Lift the BigFloat results back to `Arb` (or `Acb` for complex)
     with radius zero, so the return-type element matches input.

Same recipe for `Acb` via `Complex{BigFloat}`.

## Precision caveat

The `BigFloat` precision used for the SVD is **the current
`setprecision(BigFloat)`**, NOT each `Arb` value's per-value
precision.  For consistent results, set
`setprecision(BigFloat, p)` to match the precision of your `Arb`
inputs before calling `pade_svd`.  A future v2 could read each
input's precision and call `setprecision` per-call; v1 takes the
simpler global-precision approach.

## What's NOT included

  - `pade_svd` for `Acb` matrices over non-Float `Arb` precisions
    higher than the global `BigFloat` precision: clamped at the
    global precision (same caveat as above).
  - `Polynomials.roots(::Polynomial{Arb})` — separate friction bead
    `padetaylor-8pi`; affects `StepControl.step_pade_root` only.
    Not in scope here.
  - `default_tol(::Type{Arb})` — defer to user-supplied `tol` kwarg
    for v1.

## References

  - `docs/adr/0003-extensions-pattern.md` §"PadeTaylorArblibExt.jl".
  - `docs/adr/0002-bigfloat-svd-via-genericlinalg.md` — the
    radius-discard caveat.
  - `RESEARCH.md §5.1` — Arblib SVD landscape (none, as confirmed
    by source).
  - `src/LinAlg.jl` — the module this extension specialises.
"""
module PadeTaylorArblibExt

import PadeTaylor.LinAlg
using GenericLinearAlgebra: GenericLinearAlgebra
using Arblib:               Arb, Acb

# =============================================================================
# pade_svd(::Matrix{Arb})
# =============================================================================

"""
    pade_svd(A::AbstractMatrix{Arb}; full::Bool = false)
        -> (U::Matrix{Arb}, S::Vector{Arb}, Vt::Matrix{Arb})

Route through `BigFloat`: convert `A` to `Matrix{BigFloat}` (radius
discarded), run Jacobi SVD, lift back to `Arb` with radius zero.
"""
function LinAlg.pade_svd(A::AbstractMatrix{Arb}; full::Bool = false)
    A_bf = BigFloat.(A)
    F    = GenericLinearAlgebra.svd(A_bf; full = full)
    U    = Arb.(F.U)
    S    = Arb.(F.S)
    Vt   = Arb.(F.Vt)
    return (U, S, Vt)
end

# =============================================================================
# pade_svd(::Matrix{Acb})
# =============================================================================

"""
    pade_svd(A::AbstractMatrix{Acb}; full::Bool = false)
        -> (U::Matrix{Acb}, S::Vector{Arb}, Vt::Matrix{Acb})

Route through `Complex{BigFloat}`: convert `A` to
`Matrix{Complex{BigFloat}}` (radius discarded), run Jacobi SVD, lift
the U/Vt back to `Acb`; singular values are real, lifted to `Arb`.
"""
function LinAlg.pade_svd(A::AbstractMatrix{Acb}; full::Bool = false)
    A_cbf = Complex{BigFloat}.(A)
    F     = GenericLinearAlgebra.svd(A_cbf; full = full)
    U     = Acb.(F.U)
    S     = Arb.(F.S)
    Vt    = Acb.(F.Vt)
    return (U, S, Vt)
end

end # module PadeTaylorArblibExt
