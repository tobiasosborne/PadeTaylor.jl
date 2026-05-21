"""
    PadeTaylor.VectorBVP

Chebyshev spectral-collocation boundary-value-problem solver for
**first-order vector ODE systems** on a complex line segment

    y'(z) = f(z, y),    y : [z_a, z_b] ⊂ ℂ → ℂ^d,

subject to a **general linear two-point boundary condition**

    B_a · y(z_a) + B_b · y(z_b) = g          (d scalar equations).

This module is the vector first-order analogue of the scalar
second-order `BVP.bvp_solve`.  It exists because the v0.2 headline
target — the `P_I⁽²⁾` tritronquée solution (bead `padetaylor-r2m`,
KKG 2015 §7) — is a **separatrix**: forward IVP integration diverges
(one mode grows exponentially), so a global collocation method is
*mandatory*.  Shooting cannot pin a separatrix; collocation never
integrates forward, so it is immune to the instability.

## Why a new module (not a lift of `BVP.jl`)

The scalar `BVP.bvp_solve` is hard-wired to a *second-order scalar*
ODE `u'' = f(z,u)` with *Dirichlet* boundary conditions.  The
`P_I⁽²⁾` companion form is a first-order *4-vector* system whose
boundary data is a general linear condition, not pointwise Dirichlet.
Lifting `BVP.jl` in place would break the v0.1 scalar tests or burst
the 200-LOC cap; an additive new module is the honest answer
(ADR-0023).

## The algorithm — Chebyshev collocation + Newton + τ-method

  1. **Discretise.**  `N+1` Chebyshev extrema nodes
     `t_j = cos(jπ/N)`, `j = 0…N` (DMSUITE convention: `t_0 = +1`,
     `t_N = -1`).  The affine map `z(t) = (z_a+z_b)/2 + s·t` with
     `s = (z_b-z_a)/2` sends `t = -1 ↦ z_a` and `t = +1 ↦ z_b`.  In
     the variable `t` the system becomes `dY/dt = s·f(z(t), Y)` — the
     chain rule contributes the single scale factor `s`.

  2. **`D1`.**  The `(N+1)×(N+1)` Chebyshev first-derivative matrix
     (Trefethen `cheb.m` / Weideman-Reddy `chebdif.m`).  The builder
     is copied verbatim from `BVP._chebyshev_D1` — we do *not* reach
     into another module's private symbol.

  3. **Stacked unknown.**  `Y` is the `d`-vector at each of the `N+1`
     nodes; it is solved as one stacked vector of length `(N+1)·d`,
     node-major (`[Y_0; Y_1; …; Y_N]`).  The collocation derivative
     operator on the stack is `D1 ⊗ I_d`.

  4. **Residual.**  Collocation residual at node `j`:
     `R_j = (D1·Y)_j − s·f(z_j, Y_j)` — `(N+1)·d` equations.  The BC
     contributes `d` more: `B_a·Y_0 + B_b·Y_N − g`.  The system is
     over-determined by `d`; the **τ-method** replaces the `d`
     collocation rows at node `t_0` (the `z_a` endpoint) by the `d`
     BC rows.  (If the linear oracle below floored above machine
     precision this row choice would be wrong — it does not; see
     ADR-0023 and the `VB.1.x` test.)

  5. **Newton.**  Solve `J·ΔY = R`, update `Y ← Y − ΔY`, stop at
     `‖ΔY‖_∞ ≤ tol`.  The Jacobian's collocation rows are
     `(D1 ⊗ I_d) − s·blkdiag(∂f/∂y(z_j,Y_j))`; the `d` replaced rows
     are the constant BC block `[B_a in node-0 cols | … | B_b in
     node-N cols]`.  Step-norm convergence (not residual-norm) is the
     spectral-BVP-correct test — the discrete residual floors at
     `cond(D1⊗I)·eps`, unachievable as a Newton target (same rationale
     as scalar `BVP.bvp_solve`).

  6. **Node Jacobian — `Taylor1` autodiff.**  `∂f/∂y(z_j,Y_j)` is the
     `d×d` matrix obtained *exactly* (no finite differences, no new
     dependency): promote `Y_j` to `Vector{Taylor1{CT}}` of constant
     series; for column `i`, seed component `i` with a unit linear
     series `y[i] + Taylor1([0,1])`, evaluate `f`, and read
     coefficient `[1]` of every output component.  `d` RHS evaluations
     give the full matrix.  This rides on the element-type-generic RHS
     convention of ADR-0020.  A caller may pass `jacobian = Jf`
     (`(z,y) -> d×d matrix`) to bypass autodiff.

  7. **Dense output.**  Berrut-Trefethen 2004 Chebyshev-2 barycentric
     interpolation applied component-wise; `sol(z) -> Vector{CT}`.
     Out-of-segment `z` throws `DomainError`.

## Generic-`T` discipline

Every operation routes through Julia's promote-type machinery.
Float64 takes LAPACK paths (matrix multiply, backslash); BigFloat
routes the Newton solve through `GenericLinearAlgebra` — the same
idiom as `BVP.jl`.  `Complex{T}` endpoints propagate through `s`.

## Fail-fast contract (CLAUDE.md Rule 1)

Throws with a `Suggestion` / `detail` line on: `N < 4`; `maxiter < 1`;
`tol < 0`; dimension mismatches in `B_a`/`B_b`/`g` or the RHS output;
non-convergent Newton; a singular Jacobian; out-of-segment evaluation.

## References

  - `src/BVP.jl` — the scalar second-order solver this module mirrors;
    verbatim source of the `_chebyshev_D1` builder.
  - `references/bvp_recipe.md` — §1 nodes + affine map, §2 `D1`, §5
    barycentric evaluation.
  - `docs/adr/0023-vector-bvp-solver.md` — the design decision.
  - ADR-0020 (element-type-generic vector RHS), ADR-0022 (`P_I⁽²⁾`
    companion form), ADR-0002 (BigFloat linear algebra).
"""
module VectorBVP

using LinearAlgebra: I, kron, Diagonal
using GenericLinearAlgebra: GenericLinearAlgebra  # extends `\` for BigFloat
using TaylorSeries: Taylor1

export vector_bvp_solve, VectorBVPSolution

# =============================================================================
# Solution container
# =============================================================================

"""
    VectorBVPSolution{T <: AbstractFloat, CT <: Number}

The converged Chebyshev-collocation vector BVP solution on `[z_a, z_b]`.

Fields:
  - `z_a, z_b`               : segment endpoints (complex in general).
  - `N`                      : `N+1 = length(nodes_t)` collocation nodes.
  - `nodes_t::Vector{T}`     : Chebyshev extrema in `[-1,1]`, descending.
  - `nodes_z::Vector{CT}`    : nodes after the affine map.
  - `Y_nodes::Vector{Vector{CT}}` : converged `d`-vector at each node.
  - `residual_inf::T`        : final `‖R‖_∞` of the Newton residual.
  - `iterations::Int`        : Newton iterations taken to converge.

Call as `sol(z) -> Vector{CT}` for barycentric evaluation of the `d`
components at any `z ∈ [z_a, z_b]`.
"""
struct VectorBVPSolution{T <: AbstractFloat, CT <: Number}
    z_a          :: CT
    z_b          :: CT
    N            :: Int
    nodes_t      :: Vector{T}
    nodes_z      :: Vector{CT}
    Y_nodes      :: Vector{Vector{CT}}
    residual_inf :: T
    iterations   :: Int
end

# =============================================================================
# Public driver
# =============================================================================

"""
    vector_bvp_solve(f, z_a, z_b, B_a, B_b, g;
                     N::Integer = 20, tol = nothing, maxiter::Integer = 10,
                     initial_guess = nothing, jacobian = nothing)
        -> VectorBVPSolution

Solve the first-order vector system `y'(z) = f(z, y)` on the complex
segment `[z_a, z_b]` under the general linear boundary condition
`B_a·y(z_a) + B_b·y(z_b) = g` via Chebyshev-collocation Newton.

Arguments:
  - `f(z, y) -> Vector`   : the vector RHS; `y` a `d`-vector whose
    elements may be plain `CT` or `Taylor1{CT}` (element-type generic).
  - `B_a, B_b`            : `d×d` boundary-condition matrices.
  - `g`                   : the length-`d` boundary right-hand side;
    `d` is inferred from `length(g)`.

Kwargs:
  - `N::Integer = 20`     : `N+1` collocation nodes; must be `≥ 4`.
  - `tol = nothing`       : Newton step-norm tolerance; default
    `eps(real(T))^(3/4)`.
  - `maxiter::Integer = 10` : Newton iteration cap.
  - `initial_guess = nothing` : callable `z -> Vector` (length `d`);
    default is a zero `d`-vector at every node.
  - `jacobian = nothing`  : optional `(z,y) -> d×d matrix`; when
    omitted the `d×d` node Jacobian is taken by `Taylor1` autodiff.

Throws `ArgumentError` on out-of-range parameters or dimension
mismatches and `ErrorException` on non-convergence or a singular
Jacobian (each with a suggestion).
"""
function vector_bvp_solve(f, z_a, z_b, B_a, B_b, g;
                          N::Integer = 20,
                          tol = nothing,
                          maxiter::Integer = 10,
                          initial_guess = nothing,
                          jacobian = nothing)

    # ---- validate -----------------------------------------------------------
    N ≥ 4 || throw(ArgumentError(
        "vector_bvp_solve: N must be ≥ 4 (got $N).  Suggestion: increase to " *
        "≥ 10 — the algorithm needs a handful of interior collocation nodes."))
    maxiter ≥ 1 || throw(ArgumentError(
        "vector_bvp_solve: maxiter must be ≥ 1 (got $maxiter)."))

    d = length(g)
    size(B_a) == (d, d) || throw(ArgumentError(
        "vector_bvp_solve: B_a must be $d×$d (d = length(g) = $d), got " *
        "$(size(B_a)).  Suggestion: B_a, B_b are d×d, g is length d."))
    size(B_b) == (d, d) || throw(ArgumentError(
        "vector_bvp_solve: B_b must be $d×$d (d = length(g) = $d), got " *
        "$(size(B_b))."))

    # ---- type promotion -----------------------------------------------------
    CT = promote_type(typeof(z_a), typeof(z_b),
                       eltype(B_a), eltype(B_b), eltype(g))
    T  = float(real(CT))
    tol_T = tol === nothing ? eps(T)^(T(3) / T(4)) : T(tol)
    tol_T ≥ 0 || throw(ArgumentError(
        "vector_bvp_solve: tol must be non-negative (got $tol_T)."))

    z_a_CT, z_b_CT = CT(z_a), CT(z_b)
    Ba = Matrix{CT}(B_a); Bb = Matrix{CT}(B_b); g_CT = Vector{CT}(g)

    # ---- Chebyshev extrema nodes + affine map -------------------------------
    pi_T    = T(π)
    nodes_t = T[cos(T(j) * pi_T / T(N)) for j in 0:N]
    half_sum  = (z_a_CT + z_b_CT) / 2
    s         = (z_b_CT - z_a_CT) / 2          # affine scale; dY/dt = s·f
    nodes_z   = CT[half_sum + s * t for t in nodes_t]
    np1       = N + 1

    D1  = _chebyshev_D1(nodes_t, T, N)
    Dop = kron(D1, Matrix{CT}(I, d, d))        # (D1 ⊗ I_d), the stack operator

    # ---- initial guess (stacked, node-major) --------------------------------
    Y = if initial_guess === nothing
        zeros(CT, np1 * d)
    else
        stk = Vector{CT}(undef, np1 * d)
        for j in 1:np1
            yj = initial_guess(nodes_z[j])
            length(yj) == d || throw(ArgumentError(
                "vector_bvp_solve: initial_guess returned a length-$(length(yj)) " *
                "vector, expected d = $d."))
            stk[(j-1)*d + 1 : j*d] .= CT.(yj)
        end
        stk
    end

    # ---- τ-method: the d BC rows replace node-0's d collocation rows --------
    # DMSUITE node order is descending: t_0 = +1 ↦ z_b, t_N = -1 ↦ z_a (the
    # affine map z(t) = half_sum + s·t).  So `y(z_a)` is the node-N block and
    # `y(z_b)` is the node-0 block: B_a multiplies node-N columns, B_b
    # multiplies node-0 columns.  Getting this endpoint↔node pairing wrong
    # is the τ-method's subtle trap (ADR-0023; mutation B in the test).
    bc_rows = 1:d                              # node-0 stacked rows replaced
    nodeN   = N * d
    J_bc    = zeros(CT, d, np1 * d)
    J_bc[:, nodeN+1 : nodeN+d] .= Ba           # B_a hits y(z_a) = node N
    J_bc[:, 1:d]               .= Bb           # B_b hits y(z_b) = node 0

    # ---- Newton iteration ---------------------------------------------------
    iter         = 0
    residual_inf = T(Inf)
    step_inf     = T(Inf)
    converged    = false

    for k in 1:maxiter
        iter = k
        # Collocation residual R = (D1⊗I)·Y − s·f(z_j, Y_j), stacked.
        R = Dop * Y
        for j in 1:np1
            yj = Y[(j-1)*d + 1 : j*d]
            fj = f(nodes_z[j], yj)
            length(fj) == d || throw(ArgumentError(
                "vector_bvp_solve: f(z, y) returned a length-$(length(fj)) " *
                "vector, expected d = $d.  Suggestion: the RHS must return " *
                "a d-vector matching length(g)."))
            R[(j-1)*d + 1 : j*d] .-= s .* CT.(fj)
        end
        residual_inf = T(maximum(abs, R))

        # Jacobian: collocation rows = (D1⊗I) − s·blkdiag(∂f/∂y).
        J = copy(Dop)
        for j in 1:np1
            yj  = Y[(j-1)*d + 1 : j*d]
            Jfj = jacobian === nothing ?
                  _autodiff_jacobian(f, nodes_z[j], yj, d, CT) :
                  Matrix{CT}(jacobian(nodes_z[j], yj))
            rng = (j-1)*d + 1 : j*d
            J[rng, rng] .-= s .* Jfj
        end

        # τ-method row replacement: node-0's d rows become the BC rows.
        # y(z_a) is the node-N block, y(z_b) the node-0 block (see above).
        J[bc_rows, :] .= J_bc
        Yza = Y[nodeN+1 : nodeN+d]; Yzb = Y[1:d]
        R[bc_rows]    .= Ba * Yza .+ Bb * Yzb .- g_CT

        ΔY = _solve(J, R)
        step_inf = T(maximum(abs, ΔY))
        Y .-= ΔY

        if step_inf ≤ tol_T
            converged = true
            break
        end
    end

    if !converged
        throw(ErrorException(
            "vector_bvp_solve: Newton did not converge in $maxiter iterations.  " *
            "Final ‖ΔY‖_∞ = $step_inf > tol = $tol_T (‖R‖_∞ = $residual_inf).  " *
            "Suggestion: (a) increase N; (b) provide a closer `initial_guess`; " *
            "(c) raise `maxiter`; (d) verify the BVP has a unique solution."))
    end

    # Final residual diagnostic at the converged Y (collocation part only).
    Rf = Dop * Y
    for j in 1:np1
        yj = Y[(j-1)*d + 1 : j*d]
        Rf[(j-1)*d + 1 : j*d] .-= s .* CT.(f(nodes_z[j], yj))
    end
    Rf[bc_rows] .= Ba * Y[nodeN+1 : nodeN+d] .+ Bb * Y[1:d] .- g_CT
    residual_inf = T(maximum(abs, Rf))

    Y_nodes = [Y[(j-1)*d + 1 : j*d] for j in 1:np1]
    return VectorBVPSolution{T, CT}(z_a_CT, z_b_CT, N, nodes_t, nodes_z,
                                    Y_nodes, residual_inf, iter)
end

# =============================================================================
# Callable: barycentric vector evaluation at z* ∈ [z_a, z_b]
# =============================================================================

"""
    (sol::VectorBVPSolution)(z) -> Vector{CT}

Evaluate the converged vector BVP solution at an arbitrary `z` in the
segment `[sol.z_a, sol.z_b]` via component-wise Berrut-Trefethen 2004
Chebyshev-2 barycentric interpolation.  Throws `DomainError` for `z`
outside the segment (silent extrapolation on a spectral interpolant is
unsafe).
"""
function (sol::VectorBVPSolution{T, CT})(z) where {T, CT}
    z_CT      = CT(z)
    half_diff = (sol.z_b - sol.z_a) / 2
    t_star    = (z_CT - (sol.z_a + sol.z_b) / 2) / half_diff

    real(t_star) ≤ 1 + 100*eps(T) && real(t_star) ≥ -1 - 100*eps(T) ||
        throw(DomainError(z,
            "VectorBVPSolution: z = $z is outside [z_a, z_b] = " *
            "[$(sol.z_a), $(sol.z_b)] (t* = $t_star).  Suggestion: restrict " *
            "the query to the segment."))

    d = length(sol.Y_nodes[1])
    return CT[_barycentric_eval(sol.nodes_t,
                                CT[yk[i] for yk in sol.Y_nodes], t_star)
              for i in 1:d]
end

# =============================================================================
# Helpers (un-exported; module-private)
# =============================================================================

"""
    _chebyshev_D1(t, ::Type{T}, N) -> Matrix{T}

The `(N+1)×(N+1)` Chebyshev first-derivative matrix in `t ∈ [-1,1]`,
Trefethen `cheb.m` / Weideman-Reddy `chebdif.m` formula.  Copied
verbatim from `BVP._chebyshev_D1` per ADR-0023 (a sibling module's
private symbol must not be reached into).
"""
function _chebyshev_D1(t::AbstractVector{T}, ::Type{T}, N::Integer) where T <: AbstractFloat
    np1 = N + 1
    D   = Matrix{T}(undef, np1, np1)
    @inbounds for j in 1:np1, i in 1:np1
        if i == j
            D[i, j] = zero(T)
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
    _autodiff_jacobian(f, z, y, d, ::Type{CT}) -> Matrix{CT}

The `d×d` Jacobian `∂f/∂y(z, y)` by `Taylor1` automatic differentiation.
For column `i` the `i`-th component of `y` is promoted to a unit linear
series `y[i] + Taylor1([0,1])`, the others to constant series; `f` is
evaluated and coefficient `[1]` of each output component is read off.
Exact to machine precision; rides on the element-type-generic RHS
convention (ADR-0020).
"""
function _autodiff_jacobian(f, z, y, d::Integer, ::Type{CT}) where CT
    Jf = Matrix{CT}(undef, d, d)
    @inbounds for i in 1:d
        yt = Taylor1{CT}[Taylor1(CT[y[k]], 1) for k in 1:d]
        yt[i] = Taylor1(CT[y[i], one(CT)], 1)
        ft = f(z, yt)
        length(ft) == d || throw(ArgumentError(
            "vector_bvp_solve: f(z, y) returned a length-$(length(ft)) " *
            "vector under autodiff, expected d = $d."))
        for r in 1:d
            Jf[r, i] = ft[r][1]
        end
    end
    return Jf
end

"""
    _solve(J, R) -> ΔY

Newton linear solve `J·ΔY = R`, generic in element type.  Float64
takes the LAPACK `\\` path; BigFloat routes through `GenericLinearAlgebra`
(loaded above).  A singular `J` is re-thrown as an `ErrorException` with
a suggestion (Rule 1) — `vector_bvp_solve` never returns a NaN lie.
"""
function _solve(J, R)
    try
        return J \ R
    catch err
        throw(ErrorException(
            "vector_bvp_solve: the Newton Jacobian is singular " *
            "($(typeof(err))).  Suggestion: (a) the boundary condition " *
            "B_a·y_a + B_b·y_b = g may be rank-deficient — check B_a, B_b; " *
            "(b) the linearised problem may be ill-posed at this iterate; " *
            "(c) try a different `initial_guess`."))
    end
end

"""
    _barycentric_eval(t_nodes, vals, t_star) -> CT

Berrut-Trefethen 2004 second-form barycentric interpolation on
Chebyshev-2 nodes (weights `w_j = (-1)^j`, halved at the endpoints).
The coincident-node guard returns `vals[k]` directly when `t_star`
lands on a node.  Verbatim in shape from `BVP._barycentric_eval`.
"""
function _barycentric_eval(t_nodes::AbstractVector{T},
                           vals::AbstractVector{CT},
                           t_star) where {T <: AbstractFloat, CT <: Number}
    np1  = length(t_nodes)
    atol = eps(real(CT))
    @inbounds for j in 1:np1
        abs(t_star - t_nodes[j]) ≤ atol && return vals[j]
    end
    num = zero(CT)
    den = zero(CT)
    @inbounds for j in 1:np1
        sign_w = iseven(j - 1) ? T(1) : T(-1)
        half   = (j == 1 || j == np1) ? T(0.5) : T(1)
        w      = sign_w * half
        delta  = t_star - CT(t_nodes[j])
        num   += w * vals[j] / delta
        den   += w           / delta
    end
    return num / den
end

end # module VectorBVP
