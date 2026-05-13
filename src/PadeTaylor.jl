"""
    PadeTaylor

A general-purpose Taylor–Padé IVP solver for analytic ODEs `y' = f(z, y)`
with `z ∈ ℂ`, in the Fornberg–Weideman tradition.

## Architecture

Four algorithmically independent layers, per ADR-0001:

  1. `LinAlg`        — SVD dispatcher with relative-accuracy guarantees on
                       small singular values (per ADR-0002).
  2. `RobustPade`    — GGT 2013 Algorithm 2 + Chebfun reweighting; the
                       single point of truth for Padé conversion.
  3. `Coefficients`  — Wraps `TaylorSeries.jl::Taylor1{T}` for our use.
  4. `StepControl`   — Jorba-Zou 2005 §3.2 (default) and FW 2011 §3.1
                       Padé-root distance (alternative).
  5. `PadeStepper`   — Orchestrates one step: Taylor → Padé → step.
  6. `Problems`      — `PadeTaylorProblem` / `solve_pade` public API.
  7. `PathNetwork`   — FW 2011 §3.1 5-direction wedge path-tree (Tier-2).
  8. `BVP`           — Chebyshev-Newton spectral BVP solver (Tier-3).
  9. `Dispatcher`    — 1D IVP↔BVP chain composition per FW 2011 §4.4.
 10. `EdgeDetector`  — 5-point Laplacian pole-field classifier (FW §3.2.2).

## Determinism

Same input bytes + same package version + same `T` → bit-identical output
bytes for `T <: AbstractFloat` symbolic-tier operations (none of the
core algorithm uses floating-point reductions over hash-set orderings).

Numerical-tier operations (the SVD and Padé root-finding) are bit-
identical *given the platform fingerprint* `{arch, os, runtime}` only;
this matches the workbench's `numerical: true` determinism tier
(ADR-0015 in `scientist-workbench`).

For arb-prec runs (`T = BigFloat` or `T = Arblib.Arb`), the user-set
precision is part of the input identity; same precision → same output
bytes.

## References

  - Fornberg & Weideman, *A numerical methodology for the Painlevé
    equations*, J. Comput. Phys. 230 (2011) 5957–5973. The foundational
    paper. See `references/markdown/FW2011_*.md`.
  - Gonnet, Güttel & Trefethen, *Robust Padé Approximation via SVD*,
    SIAM Review 55 (2013) 101–117. The Padé routine. See
    `references/markdown/GGT2013_*.md`.
  - Jorba & Zou, *A software package for the numerical integration of
    ODE by means of high-order Taylor methods*, Experimental
    Mathematics 14 (2005) 99–117. Step-size formula.

Full design rationale in `DESIGN.md`; deep research in `RESEARCH.md`.
"""
module PadeTaylor

# Internal modules — these compose into the public API; not re-exported.
include("LinAlg.jl")
include("RobustPade.jl")
include("Coefficients.jl")
include("StepControl.jl")
include("PadeStepper.jl")
include("Problems.jl")
include("PathNetwork.jl")
include("BVP.jl")
include("Dispatcher.jl")
include("EdgeDetector.jl")

# Public API (re-exported from sub-modules).
using .Problems:    PadeTaylorProblem, solve_pade, PadeTaylorSolution, taylor_eval
using .RobustPade:  robust_pade, PadeApproximant
using .Coefficients: taylor_coefficients_1st, taylor_coefficients_2nd
using .PathNetwork: path_network_solve, PathNetworkSolution
using .BVP:         bvp_solve, BVPSolution
using .Dispatcher:  dispatch_solve, DispatcherSolution, IVPSegment, BVPSegment
using .EdgeDetector: laplacian_residual, pole_field_mask

# CommonSolve adapter: the algorithm struct is declared HERE in the main
# module so users can construct it after `using PadeTaylor, CommonSolve`
# without a qualified name (the ext file adds init/step!/solve! methods
# on it).  Per ADR-0003 "Translation only — no algorithmic logic in the
# extension".  The methods live in `ext/PadeTaylorCommonSolveExt.jl`.
"""
    PadeTaylorAlg{H <: Real}(; h_max, max_steps = 100_000)

`CommonSolve.jl`-compatible algorithm marker for `solve(prob, alg)` /
`init(prob, alg)`.  Loaded only when `using CommonSolve` activates the
`PadeTaylorCommonSolveExt` extension; constructing without that load
gives a plain struct with no integrator methods attached.
"""
struct PadeTaylorAlg{H <: Real}
    h_max     :: H
    max_steps :: Int
end
PadeTaylorAlg(; h_max::Real, max_steps::Integer = 100_000) =
    PadeTaylorAlg{typeof(h_max)}(h_max, Int(max_steps))

export PadeTaylorProblem, solve_pade, PadeTaylorSolution, taylor_eval
export robust_pade, PadeApproximant
export taylor_coefficients_1st, taylor_coefficients_2nd
export path_network_solve, PathNetworkSolution
export bvp_solve, BVPSolution
export dispatch_solve, DispatcherSolution, IVPSegment, BVPSegment
export laplacian_residual, pole_field_mask
export PadeTaylorAlg

end # module PadeTaylor
