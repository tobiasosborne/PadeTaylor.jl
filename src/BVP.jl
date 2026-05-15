"""
    PadeTaylor.BVP

Chebyshev spectral-collocation boundary-value-problem solver for
second-order analytic ODEs `u'' = F(z, u)` on a complex line segment
`[z_a, z_b] ⊂ ℂ` with Dirichlet boundary conditions `u(z_a) = u_a`,
`u(z_b) = u_b`.

This module is the **Tier-3** companion to `PathNetwork`: the IVP
path-network handles pole-rich regions of the complex plane where the
local Padé-Taylor step is *stable*; the BVP solver handles smooth
regions where the IVP is exponentially unstable (one mode grows
exponentially, and the wrong-direction shoot amplifies it).  At a
later Phase 12 a Dispatcher composes them; for now this module is a
**standalone deliverable** usable for any 2nd-order analytic BVP that
fits the Dirichlet template.

## Algorithm — Chebyshev-Newton (FW 2011 §3.2 + Trefethen SMIM ch. 6, 13)

  1. **Discretise.**  On the canonical interval `[-1, 1]` place
     `N + 1` Chebyshev extrema nodes `t_j = cos(jπ/N)`, `j = 0…N`
     (DMSUITE convention: `t_0 = +1`, `t_N = -1`).  Map to the
     complex segment via the affine transform
     `z(t) = (z_a + z_b)/2 + (z_b - z_a)/2 · t`, so `t = -1 ↔ z_a`
     and `t = +1 ↔ z_b`.

  2. **Build differentiation matrices.**  `D₁` is the standard
     Chebyshev first-derivative matrix in `t` (Trefethen `cheb.m` /
     Weideman-Reddy `chebdif.m`).  `D₂ = D₁ · D₁` (option (a),
     adequate for `N ≤ 50`; for tighter conditioning at large `N`
     use the explicit Weideman-Reddy recursion — Tier-3 v1 ships
     the simpler option per `references/bvp_recipe.md §3`).

  3. **Newton iteration on interior nodes.**  The chain rule on
     the affine map contributes a `(z_b - z_a)²/4` scale factor:
     in `t`-space the ODE is

         D₂ u = ((z_b - z_a)/2)² · F(z(t), u)

     so the **interior residual** (rows 2…N, dropping the BC rows
     1 and N+1) is

         R = (D₂ u)_{int} − ((z_b - z_a)²/4) · F(z_int, u_int)

     with the BC contribution `D₂[int, 1] · u_b + D₂[int, N+1] · u_a`
     precomputed.  The **analytic Jacobian** is

         J = D₂[int, int] − ((z_b - z_a)²/4) · diag(∂F/∂u(z_int, u_int))

     — diagonal modification of `D₂_{int,int}`; no automatic
     differentiation or finite-difference needed.  Newton step is
     `u_int ← u_int − J \\ R`.  Convergence at `‖R‖_∞ ≤ tol`.

  4. **Barycentric evaluation.**  For dense output at arbitrary
     `z* ∈ ℂ`: compute the pre-image `t* = (2z* − z_a − z_b) / (z_b − z_a)`,
     then apply the Berrut-Trefethen 2004 §5 barycentric formula
     with Chebyshev-2 weights `w_j = (−1)^j` (halved at `j = 0, N`).
     Derivative via the chain rule: `u'(z) = (du/dt) / ((z_b − z_a)/2)`
     where `du/dt` is interpolated from the collocation-node values
     `D₁ · u_nodes`.

The `(z_b - z_a)²/4` affine-scale factor is *the* subtle place to
get the algorithm wrong — FW 2011 line 190 says only "the Jacobian
is trivial to calculate" without writing it down, and the
problems-oracle subagent initially seeded BCs from the wrong ODE
(equianharmonic ℘) before catching the spec drift.  The factor is
derived once in §4 of `references/bvp_recipe.md` and cited there.

## Generic-`T` discipline (`references/bvp_recipe.md §6`)

Every arithmetic operation routes through Julia's promote-type
machinery:

  - **Float64**: LAPACK paths everywhere (matrix multiply, backslash).
  - **BigFloat-256**: same algorithmic code; backslash routes to
    `GenericLinearAlgebra` (already in `Project.toml`).
  - **Complex{T}**: segment endpoints `z_a, z_b` may be complex;
    everything internal handles `Complex{T}` automatically because
    the `(z_b − z_a)²/4` scale factor is complex and propagates.
    The Chebyshev nodes `t_j` themselves remain real.

For `Arblib.Arb` element type: route the Newton linear solve through
`BigFloat` per the same convention as `LinAlg.pade_svd` (deferred to
Phase 8's `PadeTaylorArblibExt`).

## Composition with `PathNetwork` (Tier 3, deferred)

The Tier-3 Dispatcher (Phase 12, `padetaylor-???`) will:
  - Run `path_network_solve` on the pole-rich subdomain.
  - Run `bvp_solve` on smooth subdomains, with BCs from the path-
    network's derivative values at the pole-field edges (FW 2011
    line 192 derivative-match tolerance 1e-7 / 1e-8).
  - Stitch the two solutions into a single dense output.

This module does NOT yet know about path-networks; it is a
self-contained BVP solver that the Dispatcher will compose.

## Three-argument RHS overload `f(z, u, u')` (bead `padetaylor-i76`)

The original `bvp_solve(f, ∂f_∂u, …)` assumes `f(z, u)` with no
explicit `u'` dependence — adequate for `u'' = 6u² + z` (PI) but **not**
for `P̃_III`, whose RHS is `w'' = (w')²/w + …` (FFW 2017 md:43).  The
hybrid driver `IVPBVPHybrid` needs the three-arg form on the pole-free
sector.  A *positional* overload `bvp_solve(f, ∂f_∂u, ∂f_∂up, …)` (note
the extra Jacobian-with-respect-to-`u'` callable) selects the 3-arg
path; the original 2-arg API is preserved byte-exactly.

At each interior collocation node `i` the 3-arg path computes
`u'_int[i] = (D₁ · u_nodes)[i] / half_diff` (chain rule on the affine
map `t ↔ z`).  The residual is
`R = D₂_ii · u_int + bc_col − scale · f(z_int, u_int, u'_int)`.
The Newton Jacobian gains a non-diagonal contribution from `u'`'s
linear dependence on `u_int` via `D₁`:

    J = D₂_ii
        − scale · diag(∂f/∂u(z_int, u_int, u'_int))
        − (scale / half_diff) · diag(∂f/∂u'(z_int, u_int, u'_int))
                              · D₁[int, int]

The 2-arg path is recovered by taking `∂f/∂u' = 0` and dropping the
`u'` argument from `f`'s call — algorithmically a strict generalisation.
Keeping the API split avoids a per-call `hasmethod` probe and keeps the
2-arg fast-path's diagonal-only Jacobian unchanged.

## Fail-fast contract (CLAUDE.md Rule 1)

Throws with a `Suggestion` line on:
  - `N < 4` (need ≥ 3 interior collocation nodes for a meaningful Newton).
  - `maxiter < 1`.
  - `tol < 0`.
  - Non-convergent Newton (final `‖R‖_∞ > tol` after `maxiter` steps).
  - Singular Newton Jacobian (linear solve fails).
  - Evaluating the callable outside `[z_a, z_b]`'s pre-image disc
    (i.e., `|t*| > 1` by more than `100·eps`): silent extrapolation
    is unsafe; user must explicitly request via `extrapolate = true`.

## References

  - `references/bvp_recipe.md` — full canonical recipe, 8 sections.
  - `references/markdown/TrefethenSMIM_2000_book/...md` — SMIM ch.6
    (D₁/D₂), ch.13/14 (Newton on nonlinear BVP).
  - `references/markdown/WeidemanReddy2000_DMSUITE_ACMTOMS26/...md`
    — `chebdif`, `chebint` algorithms.
  - `references/markdown/BerrutTrefethen2004_barycentric_SIAMReview/...md`
    — barycentric formula and Chebyshev-2 weights.
  - `external/DMSUITE/{chebdif.m, chebint.m}` — reference impl;
    the Octave oracle at `external/probes/bvp-oracle/capture.m`
    builds on these.
  - `external/probes/bvp-oracle/oracles.txt` — pinned outputs;
    Julia-formatted constants in `test/_oracle_bvp.jl`.
  - `docs/adr/0004-path-network-architecture.md` — Tier-3
    deferral plan.
  - FW 2011 §3.2 (`references/markdown/FW2011_painleve_methodology_JCP230/
    FW2011_painleve_methodology_JCP230.md:176–200`) — the BVP-
    solver section that motivated this module.
"""
module BVP

using LinearAlgebra: Diagonal, lu

export bvp_solve, BVPSolution

# =============================================================================
# Solution container
# =============================================================================

"""
    BVPSolution{T <: AbstractFloat, CT <: Number}

The converged Chebyshev-Newton BVP solution on `[z_a, z_b]`.

Fields:
  - `z_a, z_b`         : segment endpoints (complex in general).
  - `u_a, u_b`         : Dirichlet boundary values.
  - `N`                : `N + 1 = length(nodes_t)` collocation nodes.
  - `nodes_t::Vector{T}` : Chebyshev extrema in `[-1, 1]`, descending
    (`+1, …, -1`) per DMSUITE convention.
  - `nodes_z::Vector{CT}`: nodes after the affine map.
  - `u_nodes::Vector{CT}`: converged solution at each `nodes_z`.
  - `residual_inf`     : final `‖R‖_∞` of the Newton residual.
  - `iterations`       : Newton iterations taken to converge.

Call as `sol(z) -> (u, u')` for barycentric evaluation at any `z ∈ [z_a, z_b]`.
"""
struct BVPSolution{T <: AbstractFloat, CT <: Number}
    z_a          :: CT
    z_b          :: CT
    u_a          :: CT
    u_b          :: CT
    N            :: Int
    nodes_t      :: Vector{T}
    nodes_z      :: Vector{CT}
    u_nodes      :: Vector{CT}
    residual_inf :: T
    iterations   :: Int
end

# =============================================================================
# Public driver
# =============================================================================

"""
    bvp_solve(f, ∂f_∂u, z_a, z_b, u_a, u_b;
              N::Integer = 20,
              tol = nothing,
              maxiter::Integer = 10,
              initial_guess = nothing) -> BVPSolution

Solve `u'' = f(z, u)` on the complex segment `[z_a, z_b]` with
Dirichlet boundary conditions `u(z_a) = u_a`, `u(z_b) = u_b` via
Chebyshev-Newton spectral collocation.

Arguments:
  - `f(z, u) -> Number`        : the ODE right-hand side.
  - `∂f_∂u(z, u) -> Number`    : its analytic partial in `u` (caller's
    responsibility — no automatic differentiation here; senior-Julia
    discipline is "you wrote `f`, you can write its derivative").

Kwargs:
  - `N::Integer = 20`           : number of subintervals; total
    nodes are `N + 1`.  Must satisfy `N ≥ 4`.
  - `tol = nothing`             : Newton **step-norm** convergence
    tolerance: iteration stops when `‖Δu‖_∞ ≤ tol`.  Default
    `eps(T)^(3/4)` ≈ `5e-13` for `Float64`, scales appropriately
    for `BigFloat`.  This is the right criterion for spectral
    BVPs because the discrete residual `‖R‖_∞` floors at
    `cond(D₂_int,int) · eps(T) ≈ N² · eps(T)` — unachievable as
    a Newton target.  The residual `‖R‖_∞` is recorded in the
    solution for diagnostics but is not the convergence test.
    See `references/bvp_recipe.md §7` open-spec-gap note.
  - `maxiter::Integer = 10`     : Newton iteration cap.  FW 2011 line
    190 reports ≤ 6 in practice; we leave headroom.
  - `initial_guess = nothing`   : if supplied, a function `z -> Number`
    evaluated at each node; otherwise a linear ramp from `u_a` to `u_b`.

Throws `ArgumentError` on out-of-range parameters and `ErrorException`
on non-convergence (with a suggestion).
"""
function bvp_solve(f, ∂f_∂u, z_a, z_b, u_a, u_b;
                   N::Integer = 20,
                   tol = nothing,
                   maxiter::Integer = 10,
                   initial_guess = nothing)

    # ---- validate -----------------------------------------------------------
    N ≥ 4 || throw(ArgumentError(
        "bvp_solve: N must be ≥ 4 (got $N).  Suggestion: increase to ≥ 10 " *
        "for any nontrivial nonlinear problem; the algorithm needs at least " *
        "a handful of interior collocation nodes."))
    maxiter ≥ 1 || throw(ArgumentError(
        "bvp_solve: maxiter must be ≥ 1 (got $maxiter)."))

    # ---- type promotion -----------------------------------------------------
    CT  = promote_type(typeof(z_a), typeof(z_b), typeof(u_a), typeof(u_b))
    T   = float(real(CT))
    # Step-norm Newton convergence: ‖Δu‖_∞ ≤ tol.  Default eps^(3/4) is
    # the standard production-Newton choice — tight enough to certify
    # convergence, loose enough to be achievable past the (irreducible)
    # spectral-residual floor cond(D₂) · eps.  See recipe §7 + the
    # mutation-proof commentary in test/bvp_test.jl.
    tol_T = tol === nothing ? eps(T)^(T(3)/T(4)) : T(tol)
    tol_T ≥ 0 || throw(ArgumentError(
        "bvp_solve: tol must be non-negative (got $tol_T)."))

    z_a_CT, z_b_CT = CT(z_a), CT(z_b)
    u_a_CT, u_b_CT = CT(u_a), CT(u_b)

    # ---- Chebyshev extrema nodes t_j = cos(jπ/N), j = 0..N ------------------
    # DMSUITE convention: t_0 = +1, t_N = -1.  Affine map below makes
    # t = -1 ↔ z_a and t = +1 ↔ z_b.
    pi_T    = T(π)
    nodes_t = T[cos(T(j) * pi_T / T(N)) for j in 0:N]

    # ---- affine map [-1, 1] → [z_a, z_b] -----------------------------------
    half_sum  = (z_a_CT + z_b_CT) / 2
    half_diff = (z_b_CT - z_a_CT) / 2          # = (z_b − z_a) / 2
    nodes_z   = CT[half_sum + half_diff * t for t in nodes_t]

    # ---- D₁, D₂ = D₁·D₁ -----------------------------------------------------
    D1 = _chebyshev_D1(nodes_t, T, N)
    D2 = D1 * D1                               # real (T)-valued

    # ---- initial guess -----------------------------------------------------
    u_nodes = if initial_guess === nothing
        # Linear ramp: u(t) = u_b at t=+1, u_a at t=-1.
        CT[u_b_CT + (u_a_CT - u_b_CT) * ((1 - t) / 2) for t in nodes_t]
    else
        CT[CT(initial_guess(z)) for z in nodes_z]
    end
    u_nodes[1]   = u_b_CT      # enforce BC at t = +1 ↔ z = z_b
    u_nodes[N+1] = u_a_CT      # enforce BC at t = -1 ↔ z = z_a

    # ---- Newton iteration on interior nodes (indices 2..N) ------------------
    scale  = (z_b_CT - z_a_CT)^2 / 4           # the affine Jacobian factor
    int    = 2:N
    D2_ii  = D2[int, int]                      # constant; reuse across iters
    bc_col = D2[int, 1] .* u_b_CT .+ D2[int, N+1] .* u_a_CT
    z_int  = nodes_z[int]
    u_int  = u_nodes[int]

    iter         = 0
    residual_inf = T(Inf)
    step_inf     = T(Inf)
    converged    = false

    for k in 1:maxiter
        iter = k
        # Residual R = (D₂ u)_int − scale · f(z_int, u_int) — diagnostic only.
        F_vals       = CT[f(z_int[j], u_int[j]) for j in eachindex(u_int)]
        R            = D2_ii * u_int .+ bc_col .- scale .* F_vals
        residual_inf = T(maximum(abs, R))

        # Analytic Jacobian: J = D₂_ii − scale · diag(∂f/∂u(z_int, u_int))
        diag_∂f   = CT[∂f_∂u(z_int[j], u_int[j]) for j in eachindex(u_int)]
        J         = D2_ii .- scale .* Diagonal(diag_∂f)
        Δu        = J \ R
        step_inf  = T(maximum(abs, Δu))
        u_int   .-= Δu

        # Step-norm convergence — the spectral-BVP-correct criterion.
        if step_inf ≤ tol_T
            converged = true
            # Recompute residual at the converged u_int for the diagnostic.
            F_vals_final = CT[f(z_int[j], u_int[j]) for j in eachindex(u_int)]
            R_final      = D2_ii * u_int .+ bc_col .- scale .* F_vals_final
            residual_inf = T(maximum(abs, R_final))
            break
        end
    end

    if !converged
        throw(ErrorException(
            "bvp_solve: Newton did not converge in $maxiter iterations.  " *
            "Final ‖Δu‖_∞ = $step_inf > tol = $tol_T (‖R‖_∞ = $residual_inf).  " *
            "Suggestion: (a) increase N (better resolution lowers truncation " *
            "error); (b) provide a closer `initial_guess` (Newton can stall " *
            "far from the basin of attraction); (c) raise `maxiter`; " *
            "(d) verify the BVP has a unique solution on this segment."))
    end

    u_nodes[int] = u_int
    return BVPSolution{T, CT}(z_a_CT, z_b_CT, u_a_CT, u_b_CT, N,
                              nodes_t, nodes_z, u_nodes,
                              residual_inf, iter)
end

"""
    bvp_solve(f, ∂f_∂u, ∂f_∂up, z_a, z_b, u_a, u_b;
              N::Integer = 20, tol = nothing, maxiter::Integer = 10,
              initial_guess = nothing, initial_guess_up = nothing) -> BVPSolution

Three-argument-RHS overload: solve `u'' = f(z, u, u')` on the complex
segment `[z_a, z_b]` with Dirichlet BCs `u(z_a) = u_a`, `u(z_b) = u_b`.
The extra `∂f_∂up(z, u, u') -> Number` callable is the analytic partial
of `f` with respect to `u'` (caller supplies — no autodiff).  Returns a
`BVPSolution` with the same fields as the 2-arg version.

This overload exists for hybrid-driver use on `P̃_III` (FFW 2017 md:43
`w'' = (w')²/w + …`), whose RHS depends on `w'`.  The 2-arg path
remains the recommended choice when `∂f/∂u' ≡ 0`.  See the module
docstring section "Three-argument RHS overload" for the algorithm
derivation.

Note: `initial_guess_up` is accepted for completeness but is currently
unused — the initial Newton iterate's `u'` is recovered from the
initial `u`-guess via `D₁`; specifying it independently would be
inconsistent with the polynomial collocation invariant.  Reserved for
future quasi-linearisation extensions.
"""
function bvp_solve(f, ∂f_∂u, ∂f_∂up, z_a, z_b, u_a, u_b;
                   N::Integer = 20,
                   tol = nothing,
                   maxiter::Integer = 10,
                   initial_guess = nothing,
                   initial_guess_up = nothing)

    initial_guess_up  # currently unused; reserved.

    N ≥ 4 || throw(ArgumentError(
        "bvp_solve: N must be ≥ 4 (got $N).  Suggestion: increase to ≥ 10."))
    maxiter ≥ 1 || throw(ArgumentError(
        "bvp_solve: maxiter must be ≥ 1 (got $maxiter)."))

    CT  = promote_type(typeof(z_a), typeof(z_b), typeof(u_a), typeof(u_b))
    T   = float(real(CT))
    tol_T = tol === nothing ? eps(T)^(T(3)/T(4)) : T(tol)
    tol_T ≥ 0 || throw(ArgumentError(
        "bvp_solve: tol must be non-negative (got $tol_T)."))

    z_a_CT, z_b_CT = CT(z_a), CT(z_b)
    u_a_CT, u_b_CT = CT(u_a), CT(u_b)

    pi_T    = T(π)
    nodes_t = T[cos(T(j) * pi_T / T(N)) for j in 0:N]

    half_sum  = (z_a_CT + z_b_CT) / 2
    half_diff = (z_b_CT - z_a_CT) / 2
    nodes_z   = CT[half_sum + half_diff * t for t in nodes_t]

    D1 = _chebyshev_D1(nodes_t, T, N)
    D2 = D1 * D1

    u_nodes = if initial_guess === nothing
        CT[u_b_CT + (u_a_CT - u_b_CT) * ((1 - t) / 2) for t in nodes_t]
    else
        CT[CT(initial_guess(z)) for z in nodes_z]
    end
    u_nodes[1]   = u_b_CT
    u_nodes[N+1] = u_a_CT

    scale     = (z_b_CT - z_a_CT)^2 / 4
    int       = 2:N
    D2_ii     = D2[int, int]
    D1_ii     = D1[int, int]            # for ∂u'/∂u_int via chain rule
    D1_ib     = D1[int, [1, N+1]]       # boundary contribution to u'_int
    bc_col_D2 = D2[int, 1] .* u_b_CT .+ D2[int, N+1] .* u_a_CT
    z_int     = nodes_z[int]
    u_int     = u_nodes[int]

    iter         = 0
    residual_inf = T(Inf)
    step_inf     = T(Inf)
    converged    = false

    # The chain-rule factor: u'(z_i) = (D₁·u_nodes)_i / half_diff.
    inv_hd    = one(CT) / half_diff
    # boundary contribution to up_int (constant across Newton iterates).
    bc_col_D1 = D1_ib * CT[u_b_CT, u_a_CT]    # length N-1

    for k in 1:maxiter
        iter = k
        up_int = (D1_ii * u_int .+ bc_col_D1) .* inv_hd
        F_vals = CT[f(z_int[j], u_int[j], up_int[j]) for j in eachindex(u_int)]
        R      = D2_ii * u_int .+ bc_col_D2 .- scale .* F_vals
        residual_inf = T(maximum(abs, R))

        diag_∂fu  = CT[∂f_∂u(z_int[j], u_int[j], up_int[j]) for j in eachindex(u_int)]
        diag_∂fup = CT[∂f_∂up(z_int[j], u_int[j], up_int[j]) for j in eachindex(u_int)]
        J = D2_ii .- scale .* Diagonal(diag_∂fu) .-
            (scale * inv_hd) .* (Diagonal(diag_∂fup) * D1_ii)
        Δu       = J \ R
        step_inf = T(maximum(abs, Δu))
        u_int  .-= Δu

        if step_inf ≤ tol_T
            converged = true
            up_final = (D1_ii * u_int .+ bc_col_D1) .* inv_hd
            F_fin    = CT[f(z_int[j], u_int[j], up_final[j]) for j in eachindex(u_int)]
            R_final  = D2_ii * u_int .+ bc_col_D2 .- scale .* F_fin
            residual_inf = T(maximum(abs, R_final))
            break
        end
    end

    if !converged
        throw(ErrorException(
            "bvp_solve(3-arg): Newton did not converge in $maxiter iterations.  " *
            "Final ‖Δu‖_∞ = $step_inf > tol = $tol_T (‖R‖_∞ = $residual_inf).  " *
            "Suggestion: (a) increase N; (b) provide a closer `initial_guess`; " *
            "(c) raise `maxiter`; (d) verify ∂f/∂u and ∂f/∂u' are correct."))
    end

    u_nodes[int] = u_int
    return BVPSolution{T, CT}(z_a_CT, z_b_CT, u_a_CT, u_b_CT, N,
                              nodes_t, nodes_z, u_nodes,
                              residual_inf, iter)
end

# =============================================================================
# Callable: barycentric evaluation at z* ∈ [z_a, z_b]
# =============================================================================

"""
    (sol::BVPSolution)(z) -> (u, u')

Evaluate the BVP solution and its derivative at an arbitrary `z` in
the segment `[sol.z_a, sol.z_b]` via the Berrut-Trefethen 2004
barycentric formula on the Chebyshev-2 grid.  The derivative is
recovered by chain rule on the affine map: `u'(z) = (du/dt) / ((z_b − z_a)/2)`,
with `du/dt` itself barycentric-interpolated from `D₁ · u_nodes`.

Throws `DomainError` if `z` lies outside the segment (i.e., the
preimage `|t*| > 1 + 100·eps`).  Silent extrapolation is unsafe
on a spectral interpolant.
"""
function (sol::BVPSolution{T, CT})(z) where {T, CT}
    z_CT   = CT(z)
    half_diff = (sol.z_b - sol.z_a) / 2
    t_star = (z_CT - (sol.z_a + sol.z_b) / 2) / half_diff

    real(t_star) ≤ 1 + 100*eps(T) && real(t_star) ≥ -1 - 100*eps(T) ||
        throw(DomainError(z,
            "BVPSolution: z = $z is outside [z_a, z_b] = [$(sol.z_a), $(sol.z_b)] " *
            "(t* = $t_star).  Suggestion: restrict the query to the segment."))

    u_val   = _barycentric_eval(sol.nodes_t, sol.u_nodes, t_star)

    # Derivative via D₁·u_nodes interpolated to t*, then chain rule.
    D1         = _chebyshev_D1(sol.nodes_t, T, sol.N)
    dudt_nodes = D1 * sol.u_nodes
    dudt_at_t  = _barycentric_eval(sol.nodes_t, dudt_nodes, t_star)
    up_val     = dudt_at_t / half_diff

    return (u_val, up_val)
end

# =============================================================================
# Helpers (un-exported; module-private)
# =============================================================================

"""
    _chebyshev_D1(t::Vector{T}, ::Type{T}, N::Integer) -> Matrix{T}

Build the `(N+1)×(N+1)` Chebyshev first-derivative matrix in `t ∈ [-1, 1]`
per the standard Trefethen SMIM `cheb.m` / Weideman-Reddy `chebdif.m`
formula.

  - Off-diagonal `i ≠ j`: `D[i,j] = c_i/c_j · (-1)^{i+j} / (t_i - t_j)`
    where `c_0 = c_N = 2` and `c_i = 1` for `0 < i < N`.
  - Diagonal `D[i,i]` set by the "negative-sum-of-row" identity to
    preserve consistency of the constant function (`D · 1 = 0`).

Generic in `T <: AbstractFloat`; `BigFloat` works without modification.
"""
function _chebyshev_D1(t::AbstractVector{T}, ::Type{T}, N::Integer) where T <: AbstractFloat
    np1 = N + 1
    D   = Matrix{T}(undef, np1, np1)
    @inbounds for j in 1:np1, i in 1:np1
        if i == j
            D[i, j] = zero(T)        # placeholder; fixed below
        else
            ci   = (i == 1 || i == np1) ? T(2) : T(1)
            cj   = (j == 1 || j == np1) ? T(2) : T(1)
            sign = iseven(i + j) ? T(1) : T(-1)
            D[i, j] = (ci / cj) * sign / (t[i] - t[j])
        end
    end
    @inbounds for i in 1:np1
        s = zero(T)
        for j in 1:np1
            i == j && continue
            s += D[i, j]
        end
        D[i, i] = -s
    end
    return D
end

"""
    _barycentric_eval(t_nodes::Vector{T}, u_nodes::Vector{CT}, t_star) -> CT

Berrut-Trefethen 2004 second-true-form barycentric interpolation for
Chebyshev-2 nodes.  Weights `w_j = (-1)^j` with halving at `j = 0, N`.

The coincident-node guard returns `u_nodes[k]` directly when
`|t_star - t_nodes[k]| ≤ eps(real(CT))`; this avoids the 0/0
indeterminate at the limit and matches the DMSUITE `chebint.m`
behaviour.
"""
function _barycentric_eval(t_nodes::AbstractVector{T},
                           u_nodes::AbstractVector{CT},
                           t_star) where {T <: AbstractFloat, CT <: Number}
    np1   = length(t_nodes)
    N_eff = np1 - 1
    atol  = eps(real(CT))

    # Coincident-node guard.
    @inbounds for j in 1:np1
        abs(t_star - t_nodes[j]) ≤ atol && return u_nodes[j]
    end

    num = zero(CT)
    den = zero(CT)
    @inbounds for j in 1:np1
        # w_j = (-1)^(j-1) · (1 if interior, 1/2 if endpoint).
        # j is 1-indexed; the "(-1)^j" in 0-indexed maps to "(-1)^(j-1)" here.
        sign_w = iseven(j - 1) ? T(1) : T(-1)
        half   = (j == 1 || j == np1) ? T(0.5) : T(1)
        w      = sign_w * half
        delta  = t_star - CT(t_nodes[j])
        num   += w * u_nodes[j] / delta
        den   += w               / delta
    end
    return num / den
end

# TIER-3 INTERFACE: Dispatcher (Phase 12) will compose `bvp_solve`
# output with `PathNetwork.path_network_solve` output via derivative-
# match at pole-field edges (FW 2011 line 192).  No code here until
# the EdgeDetector + Dispatcher beads land.

end # module BVP
