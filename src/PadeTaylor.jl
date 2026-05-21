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
  8. `Diagnostics`   — `DiagnosticReport` / `EdgeReport` containers + the
                       `quality_diagnose` generic (Delaunay-backed method
                       in `ext/PadeTaylorDiagnosticsExt.jl`); loop-closure
                       quality certificate for `PathNetworkSolution`
                       (ADR-0016).
  9. `PoleField`     — pole locations from a solved path-network's Padé
                       store (FW 2011 Fig 4.7/4.8 capability).
 10. `BVP`           — Chebyshev-Newton spectral BVP solver (Tier-3).
 11. `Dispatcher`    — 1D IVP↔BVP chain composition per FW 2011 §4.4.
 12. `EdgeDetector`  — 5-point Laplacian pole-field classifier (FW §3.2.2).
 13. `EdgeGatedSolve`    — region-growing edge-gated path-network: confines
                          the IVP to the pole field (FW §3.2.2 + md:401).
 14. `LatticeDispatcher` — 2D-lattice composition with per-row BVP fill
                          (FW §4.4); IVP source is `EdgeGatedSolve` by default
                          (ADR-0017, bead `padetaylor-0tj`).
 15. `CoordTransforms`   — Exponential coordinate maps for PIII / PV (FFW 2017 §2.1, Tier-4).
 16. `SheetTracker`      — PVI ζ-plane RHS + winding-number primitives (FFW 2017 §2.2, Tier-5).
 17. `BranchTracker`     — Walker-side cut-crossing predicate + per-branch sheet bookkeeping (ADR-0013).

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
# SharedPade lifts Padé from scalar to vector (type-II Hermite–Padé);
# it reuses LinAlg.pade_svd and RobustPade.default_tol, so it loads after
# both.  Additive — no v0.1 module behaviour changes (bead padetaylor-0ln.2).
include("SharedPade.jl")
include("Coefficients.jl")
# VectorCoefficients lifts the Taylor-jet layer from scalar to vector
# ODEs (y' = f(z, y), y ∈ ℂ^d); a d=1 jet reduces bit-identically to
# Coefficients.taylor_coefficients_1st.  Additive — no v0.1 module
# behaviour changes (bead padetaylor-0ln.7, ADR-0020).
include("VectorCoefficients.jl")
include("StepControl.jl")
include("PadeStepper.jl")
# VectorStepper takes one Padé–Taylor step of a first-order vector ODE
# y' = f(z, y), y ∈ ℂ^d, via a shared-denominator Padé (all d components
# over one Q).  It composes VectorCoefficients (vector Taylor jet) with
# SharedPade (shared-Q Hermite–Padé), so it loads after both.  Additive
# — no v0.1 module behaviour changes (bead padetaylor-0ln.8, V3a).
include("VectorStepper.jl")
# VectorStepControl is the norm-based Jorba–Zou step selector for the
# vector pipeline: it generalises StepControl.step_jorba_zou from one
# scalar jet to d jets, replacing |c_k| with a vector norm of the k-th
# coefficient vector.  At d=1 it reduces bit-identically to the scalar
# selector.  It depends only on LinearAlgebra, and VectorProblems uses
# it for the opt-in :jorba_zou adaptive policy, so it loads before
# VectorProblems.  Additive — no v0.1 module behaviour changes (bead
# padetaylor-0ln.10, V3c, ADR-0021).
include("VectorStepControl.jl")
# VectorBVP is the Chebyshev spectral-collocation Newton boundary-value
# solver for first-order vector ODE systems y' = f(z, y), y ∈ ℂ^d, under
# a general linear two-point boundary condition B_a·y(z_a) + B_b·y(z_b) =
# g.  It depends only on the standard external deps (TaylorSeries,
# LinearAlgebra, GenericLinearAlgebra) — on no sibling src/ module — so it
# could load early; it is placed BEFORE PainleveHierarchy because that
# module adds a PainleveHierarchyProblem method of vector_bvp_solve.
# Additive — the global-collocation substrate the v0.2 P_I^(2)
# tritronquée figure stands on (bead padetaylor-0ln.26, VB3, ADR-0023).
include("VectorBVP.jl")
include("Problems.jl")
# VectorProblems is the top-level driver for first-order vector ODEs
# y' = f(z, y), y ∈ ℂ^d: a problem container, a fixed-step loop over
# VectorStepper, a per-segment shared-Q store, and a dense-output
# callable.  It composes VectorStepper, so it loads after it.  Additive
# — distinct from the v0.1 solve_pade; no v0.1 dispatch changes (bead
# padetaylor-0ln.9, V3b).
include("VectorProblems.jl")
# NoumiYamada is the problem-builder layer for the even-parity
# Noumi–Yamada A_{2n}^(1) higher-order Painlevé systems: the A_{2n}^(1)
# RHS factory and the NoumiYamadaProblem builder that assembles a
# VectorPadeTaylorProblem.  It composes VectorProblems, so it loads
# after it.  Additive — no v0.1 module behaviour changes (bead
# padetaylor-0ln.12, V5a, ADR-0022).
include("NoumiYamada.jl")
# NoumiYamadaSymmetry is the affine-Weyl-group W(A_{2n}^(1)) symmetry layer
# for the Noumi–Yamada systems: the exact rational-solution oracles
# (noumi_yamada_rational) and the Bäcklund transformations
# (noumi_yamada_backlund = reflection s_i, noumi_yamada_rotation = diagram
# rotation π).  It uses NoumiYamada's RHS for its solution-preservation
# self-test, so it loads after NoumiYamada.  Additive — no v0.1/v0.2
# module behaviour changes (bead padetaylor-0ln.14, V5c).
include("NoumiYamadaSymmetry.jl")
# PainleveHierarchy is the problem-builder layer for the Painlevé-I
# hierarchy: the painleve_hierarchy(:I, m) companion-form RHS factory
# (m ∈ {1,2}; m=2 ⇒ the P_I^(2) 4-vector system), the
# PainleveHierarchyProblem builder, and the pI2_tritronquee_ic seed.
# It composes VectorProblems, so it loads after it.  Additive — no
# v0.1/v0.2 module behaviour changes (bead padetaylor-0ln.15, V6,
# ADR-0022).
include("PainleveHierarchy.jl")
# VectorPathNetworkStage2 is the Stage-2 fine-grid fill of the vector
# path-network: dense evaluation of the per-node shared-Q approximants
# over a fine grid, the substrate for the P_I^(2) tritronquée surface
# figure.  It is split out of VectorPathNetwork purely for the CLAUDE.md
# Rule 6 200-LOC cap; it has NO reverse dependency on VectorPathNetwork
# (the nearest-visited scan it needs is passed as a function argument),
# so it loads BEFORE VectorPathNetwork, which then `using`s it (bead
# padetaylor-0ln.20, F2, ADR-0015).
include("VectorPathNetworkStage2.jl")
# VectorWedgeStep is the principled wedge-step machinery of the vector
# path-network walk: the shared-Q-root-distance direction selector
# (:max_q_root) and the pole-capped adaptive step-size controller.  It
# composes StepControl (the FW 2011 §3.1 step_pade_root pole-distance
# heuristic, reused unchanged on the shared Q per ADR-0021) and
# VectorStepper, so it loads after both; VectorPathNetwork `using`s it,
# so it loads before that module.  It is split out of VectorPathNetwork
# purely for the CLAUDE.md Rule 6 200-LOC cap.  Additive — the v0.2
# B2 dense-wedge walk (bead padetaylor-0ln.37.6, absorbs 0ln.23,
# ADR-0025 Lever 2).
include("VectorWedgeStep.jl")
# VectorPathNetwork is the minimal Stage-1 vector path-network walk: a
# ‖y‖-steered 5-direction wedge that builds a visited tree of shared-Q
# approximants over a region of the complex plane.  It composes
# VectorProblems + VectorStepper, so it loads after both; it also
# `using`s VectorPathNetworkStage2 for the Stage-2 fine-grid fill, so
# that module loads before this one.  Additive — the substrate the v0.2
# A_4^(1)/P_I^(2) pole-field figures stand on (bead padetaylor-0ln.16,
# V7; Stage-2 fill bead padetaylor-0ln.20).
include("VectorPathNetwork.jl")
# VectorPoleField extracts the system's poles from each visited node's
# shared denominator Q (one Q per node — every component's poles are
# its roots).  It consumes VectorPathNetworkSolution, so it loads after
# VectorPathNetwork.  Additive — the shared-Q analogue of v0.1's
# PoleField (bead padetaylor-0ln.16, V7, ADR-0019).
include("VectorPoleField.jl")
# SheetTracker is loaded before PathNetwork because BranchTracker (the
# walker-side cut-respecting layer per ADR-0013) reuses SheetTracker's
# `winding_delta` primitive, and PathNetwork in turn uses BranchTracker.
include("SheetTracker.jl")
include("BranchTracker.jl")
# Diagnostics declares the core DiagnosticReport / EdgeReport containers
# and the empty `quality_diagnose` generic; PathNetwork then embeds the
# Union{Nothing, DiagnosticReport} field on PathNetworkSolution and
# attaches reports when path_network_solve is called with diagnose=true.
# The Delaunay-backed method lives in ext/PadeTaylorDiagnosticsExt.jl
# (ADR-0016 + ADR-0003).
include("Diagnostics.jl")
include("PathNetwork.jl")
include("PoleField.jl")
include("BVP.jl")
include("Dispatcher.jl")
include("EdgeDetector.jl")
include("EdgeGatedSolve.jl")
include("LatticeDispatcher.jl")
include("CoordTransforms.jl")
# Laplace2D is the in-house 2D-Chebyshev spectral Dirichlet Laplace
# solver on a rectangle: tensor-product Chebyshev Laplacian
# L = kron(I,D²x) + kron(D²y,I) on interior nodes, Dirichlet boundary
# data absorbed into the RHS, one linear solve.  It is voter (2) of the
# triple-method majority-vote harmonic-extension fill for the P_I^(2)
# tritronquée surface figure (voter (1): ray-fan vector_bvp_solve;
# voter (3): a Gridap.jl FEM solve homed in the figures/ project, bead
# 0ln.36).  It depends only
# on the standard external deps (LinearAlgebra, GenericLinearAlgebra) —
# on no sibling src/ module (the _chebyshev_D1 builder is copied verbatim
# from BVP.jl, not reached into).  Additive — no other module changes
# (bead padetaylor-0ln.35, ADR-0024).
include("Laplace2D.jl")
include("Painleve.jl")
include("IVPBVPHybrid.jl")
include("Heun.jl")

# Public API (re-exported from sub-modules).
using .Problems:    PadeTaylorProblem, solve_pade, PadeTaylorSolution, taylor_eval
using .RobustPade:  robust_pade, PadeApproximant
using .SharedPade:  shared_denominator_pade
using .Coefficients: taylor_coefficients_1st, taylor_coefficients_2nd
using .VectorCoefficients: vector_taylor_coefficients
using .VectorStepper: VectorPadeStepperState, vector_pade_step!,
                      vector_pade_step_with_pade!
using .VectorStepControl: vector_step_jorba_zou
using .VectorProblems: VectorPadeTaylorProblem, VectorPadeTaylorSolution,
                       vector_solve_pade
# VectorBVP — the first-order vector-system Chebyshev-collocation BVP
# solver; the global-collocation counterpart to the IVP vector_solve_pade
# (ADR-0023).
using .VectorBVP:   vector_bvp_solve, VectorBVPSolution
using .NoumiYamada: noumi_yamada_rhs, NoumiYamadaProblem
using .NoumiYamadaSymmetry: noumi_yamada_rational, noumi_yamada_backlund,
                            noumi_yamada_rotation
using .PainleveHierarchy: painleve_hierarchy, painleve_hierarchy_jacobian,
                          PainleveHierarchyProblem, pI2_tritronquee_ic
using .VectorPathNetwork: vector_path_network_solve, VectorPathNetworkSolution
using .VectorPoleField:   extract_poles_shared_q
using .PathNetwork: path_network_solve, PathNetworkSolution, eval_at, eval_at_sheet
using .Diagnostics: DiagnosticReport, EdgeReport, quality_diagnose
using .PoleField:   extract_poles
using .BVP:         bvp_solve, BVPSolution
using .Dispatcher:  dispatch_solve, DispatcherSolution, IVPSegment, BVPSegment
# Laplace2D — the in-house 2D-Chebyshev spectral Dirichlet Laplace
# solver; voter (2) of the triple-method tritronquée sector fill
# (ADR-0024).
using .Laplace2D:  laplace2d_solve, Laplace2DSolution
using .EdgeDetector: laplacian_residual, pole_field_mask
using .EdgeGatedSolve: edge_gated_pole_field_solve, EdgeGatedSolution
using .LatticeDispatcher: lattice_dispatch_solve, LatticeSolution
using .CoordTransforms: pIII_transformed_rhs, pV_transformed_rhs,
                        pIII_z_to_ζ, pIII_ζ_to_z, pV_z_to_ζ, pV_ζ_to_z
using .SheetTracker:    pVI_transformed_rhs,
                        pVI_eta_transformed_rhs, pVI_z_to_η, pVI_η_to_z,
                        winding_delta, accumulate_winding, sheet_index
using .Painleve:        PainleveProblem, PainleveSolution,
                        poles, grid_values, equation, parameters, solutionname,
                        tritronquee, hastings_mcleod
using .IVPBVPHybrid:    solve_pole_free_hybrid, IVPBVPSolution, pIII_asymptotic_ic
using .Heun:            heun_g, heun_c

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

"""
    painleveplot(sol::PainleveSolution; kwargs...) -> Makie.Figure

Plot a `PainleveSolution` in the complex `z`-plane: the solution's
trajectory (for a `solve_pade` result) or evaluation grid (for a
`path_network_solve` result), with its extracted poles overlaid.

The method is provided by the `PadeTaylorMakieExt` package extension —
it exists only when `Makie` (e.g. via `CairoMakie` / `GLMakie`) is
loaded alongside `PadeTaylor`.  Calling it without a Makie load gives a
`MethodError`.  Per ADR-0003 the extension is presentation-only: it
reads `PainleveSolution` fields and calls `poles` / `grid_values`,
nothing more.
"""
function painleveplot end

export PadeTaylorProblem, solve_pade, PadeTaylorSolution, taylor_eval
export robust_pade, PadeApproximant
export shared_denominator_pade
export taylor_coefficients_1st, taylor_coefficients_2nd
export vector_taylor_coefficients
export VectorPadeStepperState, vector_pade_step!, vector_pade_step_with_pade!
export vector_step_jorba_zou
export VectorPadeTaylorProblem, VectorPadeTaylorSolution, vector_solve_pade
export vector_bvp_solve, VectorBVPSolution
export noumi_yamada_rhs, NoumiYamadaProblem
export noumi_yamada_rational, noumi_yamada_backlund, noumi_yamada_rotation
export painleve_hierarchy, painleve_hierarchy_jacobian,
       PainleveHierarchyProblem, pI2_tritronquee_ic
export vector_path_network_solve, VectorPathNetworkSolution
export extract_poles_shared_q
export path_network_solve, PathNetworkSolution, eval_at, eval_at_sheet
export DiagnosticReport, EdgeReport, quality_diagnose
export extract_poles
export bvp_solve, BVPSolution
export dispatch_solve, DispatcherSolution, IVPSegment, BVPSegment
export laplacian_residual, pole_field_mask
export laplace2d_solve, Laplace2DSolution
export lattice_dispatch_solve, LatticeSolution
export edge_gated_pole_field_solve, EdgeGatedSolution
export pIII_transformed_rhs, pV_transformed_rhs,
       pIII_z_to_ζ, pIII_ζ_to_z, pV_z_to_ζ, pV_ζ_to_z
export pVI_transformed_rhs,
       pVI_eta_transformed_rhs, pVI_z_to_η, pVI_η_to_z,
       winding_delta, accumulate_winding, sheet_index
export PainleveProblem, PainleveSolution
export poles, grid_values, equation, parameters, solutionname
export tritronquee, hastings_mcleod
export solve_pole_free_hybrid, IVPBVPSolution, pIII_asymptotic_ic
export heun_g, heun_c
export PadeTaylorAlg
export painleveplot

end # module PadeTaylor
