"""
    PadeTaylor.SharedPadeDefect

The **relative-ODE-defect** discriminator — the per-step accuracy proxy by which
the ADR-0028 dispatcher chooses between cell A and cell B.

## Why the defect, and why it is uniform

Over one step the candidate solution in the rescaled coordinate `t ∈ [0,1]`
(`z = z₀ + h·t`, `ỹ(t) = y(z₀+h·t)`) is `ỹ_i(t) = P_i(t)/Q(t)`.  The true
solution satisfies `ỹ'(t) = h·f(z₀+h·t, ỹ(t))` (chain rule).  The **relative
defect**

    score = max_t  ‖ ỹ'(t) − h·f(z₀+h·t, ỹ(t)) ‖_∞ / (‖h·f‖_∞ + atol)

measures how badly a candidate fails its own ODE — a *reference-free* accuracy
proxy (Enright; Higham 1991; Corless–Kaya 2025, arXiv:2510.20117; the repo's own
`kkg_ode_residual` / `painleve_hierarchy_test.jl:154-186` use the same pattern).
Crucially it needs **no surrogate truth**, so — unlike the held-out-point
discriminator of ADR-0028 A1.4 — it stays valid *across* a pole: the `‖h·f‖`
normalisation cancels the leading blow-up.  The audition (probe `pole_crossing.jl`)
scored it **18/18** across pole-free / pole-crossing / BigFloat / BigFloat-crossing
— a single uniform selector with no pole detection or fallback (ADR-0028
Amendment 2).

## Two production hazards the prototype did not face

  - **RHS overflow near a pole (RISK 2).** The defect evaluates the *user* RHS at
    `y = P/Q`.  A meromorphic RHS (e.g. Calogero–Moser `1/(y₁−y₂)³`) blows up where
    `P₁=P₂` — *not* where `Q=0` — so the `Q`-root sample guard does not catch it.
    We therefore **drop any sample** whose slope or scale is non-finite; a cell that
    keeps no finite sample scores `Inf`, and the dispatcher treats two-non-finite as
    "use cell A".

  - **Extended-precision root-finding (RISK 1).** The sample guard only needs
    approximate root *locations* of `Q` to position the `t`-band.  Production runs at
    `Complex{BigFloat}`, where `Polynomials.roots` is unvalidated
    (`StepControl.jl:105-114`); so we estimate the roots in `ComplexF64` via a
    companion matrix — a deliberate precision split (the *defect arithmetic* stays at
    full `T`; only the guard-band placement is ComplexF64).  `ComplexF64` is ample to
    answer "is this real sample within `poleguard` of a pole."

## References

  - `docs/adr/0028-sharedpade-dual-construction-pareto-dispatch.md` Amendment 2 (A2.1).
  - `external/probes/adr0028-dual-construction-audition/{cells.jl,pole_crossing.jl}`
    — `relative_defect`, the 18/18 validation this module ports.
  - `src/StepControl.jl:105-114` — the `Polynomials.roots` element-type caveat behind
    the ComplexF64 root-guard split.
"""
module SharedPadeDefect

using LinearAlgebra: norm, eigvals

export relative_defect, guard_root_estimates, pole_disagreement

# Horner evaluation of a low-to-high coefficient vector at scalar `t`.
function _horner(c::AbstractVector, t)
    s = zero(t) * zero(eltype(c))
    @inbounds for k = length(c):-1:1
        s = s * t + c[k]
    end
    return s
end

# Coefficients of d/dt of a low-to-high polynomial coefficient vector.
polyder(c::AbstractVector{T}) where {T} =
    length(c) < 2 ? T[zero(T)] : T[(k - 1) * c[k] for k = 2:length(c)]

"""
    guard_root_estimates(Q) -> Vector{ComplexF64}

`ComplexF64` companion-matrix root estimates of the low-to-high coefficient vector
`Q`, used **only** to place the defect sample guard band (RISK 1: avoids
`Polynomials.roots` at `Complex{BigFloat}`).  Trailing near-zeros stripped.
"""
function guard_root_estimates(Q::AbstractVector)
    cc = ComplexF64.(Q)
    while length(cc) > 1 && abs(cc[end]) < 1e-30
        pop!(cc)
    end
    n = length(cc) - 1
    n < 1 && return ComplexF64[]
    C = zeros(ComplexF64, n, n)
    for i = 1:n-1; C[i+1, i] = 1; end
    for i = 1:n; C[i, n] = -cc[i] / cc[end]; end
    return eigvals(C)
end

"""
    relative_defect(numerators, denominator, f, z₀, h; samples, poleguard, atol)
        -> Float64

The ADR-0028 Amendment-2 selector score for one cell: the worst relative ODE
defect over interior samples `t ∈ (0,1)`, skipping a `poleguard` band around the
roots of `Q` and dropping any sample whose RHS evaluation is non-finite.  `f(z,y)`
is the ODE RHS (callable on plain numeric `y`).  Returns `Inf` if no finite sample
survives (⇒ the dispatcher uses cell A).
"""
function relative_defect(numerators::AbstractVector{<:AbstractVector{C}},
                         denominator::AbstractVector{C}, f, z₀, h;
                         samples = (0.08, 0.22, 0.37, 0.55, 0.68, 0.82, 0.94),
                         poleguard::Real = 0.05,
                         atol::Real = eps(real(C))) where {C}
    Rt = real(C)
    Q = denominator
    Qp = polyder(Q)
    numder = [polyder(p) for p in numerators]
    rts = guard_root_estimates(Q)
    worst = 0.0
    used = 0
    for tr in samples
        t = Rt(tr)
        any(r -> abs(t - r) < poleguard, rts) && continue
        qv = _horner(Q, t)
        (isfinite(qv) && abs(qv) > atol) || continue
        y = [_horner(numerators[i], t) / qv for i = 1:length(numerators)]
        qpv = _horner(Qp, t)
        yp = [(_horner(numder[i], t) * qv - _horner(numerators[i], t) * qpv) / qv^2
              for i = 1:length(numerators)]
        slope = h .* f(z₀ + h * t, y)
        all(isfinite, slope) || continue                  # RISK 2: drop overflow sample
        defect = yp .- slope
        scale = maximum(abs, slope) + atol
        (isfinite(scale) && all(isfinite, defect)) || continue
        val = Float64(maximum(abs, defect) / scale)
        isfinite(val) || continue
        worst = max(worst, val)
        used += 1
    end
    return used == 0 ? Inf : worst
end

"""
    pole_disagreement(QA, QB) -> Float64

Diagnostic: the worst `ComplexF64`-estimated root-set mismatch between two
denominators, normalised by `|root|` (ADR-0028 §5).  A large value flags an
ambiguous step.  `NaN` if either has no roots.
"""
function pole_disagreement(QA, QB)
    ra = guard_root_estimates(QA)
    rb = guard_root_estimates(QB)
    (isempty(ra) || isempty(rb)) && return NaN
    return maximum(minimum(abs(z - w) for w in rb) / max(abs(z), 1e-30) for z in ra)
end

end # module SharedPadeDefect
