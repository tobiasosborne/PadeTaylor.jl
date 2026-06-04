"""
    PadeTaylor.SharedPadeDispatch

The ADR-0028 **dual-construction dispatcher**: per vector ODE step, build both
shared-denominator Padé cells and select the more accurate one by the relative
ODE defect.  This is the layer the vector stepper calls in place of a bare
`shared_denominator_pade`, and the only place the per-step A-vs-B choice is made.

## The dispatch, precisely

`shared_pade_select(jets, m, f, z₀, h)` returns the same
`(numerators, denominator)` 2-tuple as `shared_denominator_pade`, so
`VectorStepper` destructures it unchanged.  Logic:

  1. **`d = 1` → cell A, unconditionally.** The scalar path is the bit-identity
     oracle (`SharedPade` SP.1.1, ~49 downstream tolerances); never dispatch it.
  2. **Build cell A always** (`SharedPade.shared_denominator_pade`) — the proven,
     meromorphic-safe, oracle-validated cell.
  3. **Try cell B** (`SharedPadeCellB.build_square_cell`); on any degenerate input
     it returns `nothing` ⇒ **use cell A** (`:diagonal_fallback`).
  4. **Score both** with `SharedPadeDefect.relative_defect` (needs `f`, `z₀`, `h`).
  5. **Pick** with an **A-default, type-scaled tie-break**:
     `B iff defect_B finite and defect_B < defect_A·(1 − 100·eps(real(T)))`, else A.
     A wins every tie, near-tie, and the both-non-finite case.

## Why A-default, why this tie-break (determinism, ADR C7)

The defect is a deterministic function of the jets + RHS, but a bare `≤` compare
would let an `eps`-perturbation flip the cell on a genuine close call — the exact
instability ADR-0028 A1.4 flagged against the old ε-band.  The
`(1 − 100·eps(real(T)))` band is *type-scaled* (stable at `Float64` and
`BigFloat`) and is taken **relative to the defect ratio**, not to a residual
floor.  Because a close-call misrank "costs ≈ nothing" (both cells ≈ equally
accurate, ADR-0028 Amendment 2 caveat), defaulting to A is free insurance: it
preserves the d=1, meromorphic-exactness, and FW-Table-5.1 BigFloat-A-win
contracts whenever B is not *strictly, robustly* better.

## Diagnostics (opt-in, ADR-0028 §5 / A2.3)

`diagnostics = true` returns a 3-tuple `(numerators, denominator, ::SharedPadeDispatchInfo)`
carrying `chosen_cell`, the two defect scores, the A/B pole disagreement, and the
final denominator degree.  Default `false` preserves the byte-for-byte 2-tuple
contract.

## References

  - `docs/adr/0028-sharedpade-dual-construction-pareto-dispatch.md` (Amendments 1 & 2).
  - `src/SharedPade.jl` (cell A), `src/SharedPadeCellB.jl` (cell B),
    `src/SharedPadeDefect.jl` (the selector score).
  - `src/VectorStepper.jl` — the single production call site.
"""
module SharedPadeDispatch

using ..SharedPade:       shared_denominator_pade
using ..SharedPadeCellB:  build_square_cell
using ..SharedPadeDefect: relative_defect, pole_disagreement
using ..RobustPade:       default_tol

export shared_pade_select, SharedPadeDispatchInfo

"""
    SharedPadeDispatchInfo

Per-step dispatch diagnostics (returned only when `diagnostics = true`).
`chosen_cell ∈ {:diagonal, :square, :diagonal_d1, :diagonal_fallback}`.
"""
struct SharedPadeDispatchInfo
    chosen_cell       :: Symbol
    defect_A          :: Float64
    defect_B          :: Float64
    pole_disagreement :: Float64
    denominator_degree:: Int
end

"""
    shared_pade_select(jets, m, f, z₀, h; tol, diagnostics = false)
        -> (numerators, denominator) | (numerators, denominator, ::SharedPadeDispatchInfo)

Build cell A and (for `d ≥ 2`) cell B from the rescaled `jets`, select by the
relative ODE defect (`f` = RHS, `z₀` = current point, `h` = step), and return the
chosen approximant.  See the module docstring for the d=1 short-circuit, the
fallback, and the A-default tie-break.
"""
function shared_pade_select(jets::AbstractVector{<:AbstractVector{T}}, m::Integer,
                            f, z₀, h;
                            tol::Real = default_tol(T),
                            diagnostics::Bool = false) where {T}
    d = length(jets)
    numsA, denA = shared_denominator_pade(jets, m; tol = tol)

    # d=1 short-circuit (unconditional; the scalar bit-identity contract).
    if d == 1
        return diagnostics ?
            (numsA, denA, SharedPadeDispatchInfo(:diagonal_d1, 0.0, Inf, NaN, length(denA) - 1)) :
            (numsA, denA)
    end

    cellB = build_square_cell(jets, m; tol = tol)
    if cellB === nothing
        return diagnostics ?
            (numsA, denA, SharedPadeDispatchInfo(:diagonal_fallback, NaN, NaN, NaN, length(denA) - 1)) :
            (numsA, denA)
    end
    numsB, denB = cellB

    dA = relative_defect(numsA, denA, f, z₀, h)
    dB = relative_defect(numsB, denB, f, z₀, h)

    εr = 100 * eps(real(T))
    bothbad = !isfinite(dA) && !isfinite(dB)
    bothbad && @warn "shared_pade_select: both cells scored non-finite defect; using cell A" maxlog = 5
    pickB = !bothbad && isfinite(dB) && dB < dA * (1 - εr)

    chosen = pickB ? :square : :diagonal
    nums, den = pickB ? (numsB, denB) : (numsA, denA)

    if diagnostics
        pd = pole_disagreement(denA, denB)
        return (nums, den, SharedPadeDispatchInfo(chosen, dA, dB, pd, length(den) - 1))
    end
    return (nums, den)
end

end # module SharedPadeDispatch
