# smoke.jl — validate cells.jl against known oracles before the full audition.
# Run from repo root:  julia --project=. external/probes/adr0028-dual-construction-audition/smoke.jl
#
# Checks (Rule 4 fidelity gate):
#   1. cell_A == shipped shared_denominator_pade (bit-equal mod sign/scale).
#   2. cell_A([jet],m) == robust_pade(jet,m,m;:svd);
#      cell_B_ceil([jet],m) == robust_pade(jet,m-1,m;:svd)  (d=1 oracles).
#   3. On a known rational with shared poles {2,3}, all three cells recover Q
#      and evaluate to the true f off the poles.
#   4. held_out_residual + select run and pick the exact cell on the rational.

include(joinpath(@__DIR__, "cells.jl"))
using Printf

polyval(c, z) = isempty(c) ? zero(z) : sum(c[k] * z^(k - 1) for k = 1:length(c))

# formal P/Q Taylor jet to order N
function ratio_jet(P, Q, N)
    c = zeros(ComplexF64, N + 1)
    for k = 0:N
        s = (k + 1) ≤ length(P) ? ComplexF64(P[k + 1]) : 0.0 + 0im
        for j = 1:k
            qj = (j + 1) ≤ length(Q) ? ComplexF64(Q[j + 1]) : 0.0 + 0im
            s -= qj * c[k - j + 1]
        end
        c[k + 1] = s / ComplexF64(Q[1])
    end
    c
end

function qdiff(x, y)             # ∞-norm of zero-padded difference
    L = max(length(x), length(y))
    xp = vcat(x, zeros(eltype(x), L - length(x)))
    yp = vcat(y, zeros(eltype(y), L - length(y)))
    maximum(abs.(xp .- yp))
end

println("="^70)
println("SMOKE TEST — ADR-0028 cells.jl")
println("="^70)

# --- Check 1 + 3: rational system, shared poles {2,3} -----------------------
Qt = [1.0, -5 / 6, 1 / 6]            # (1-z/2)(1-z/3), roots 2,3
P1 = [1.0, 2.0]; P2 = [3.0, -1.0]
m = 2
j1 = ratio_jet(P1, Qt, 2m + 2); j2 = ratio_jet(P2, Qt, 2m + 2)
jets = [j1, j2]

cA = cell_A(jets, m)
cBc = cell_B_ceil(jets, m)
cBf = cell_B_floor(jets, m)
nums_sh, Q_sh = shared_denominator_pade(jets, m)

# normalise Q(0)=1 already done in all; compare A to shipped
relA = qdiff(cA.denominator, Q_sh) / maximum(abs.(Q_sh))
@printf("1. cell_A vs shipped shared_denominator_pade:  Δrel(Q) = %.2e   %s\n",
        relA, relA < 1e-10 ? "OK" : "*** MISMATCH ***")
relAn = qdiff(cA.numerators[1], nums_sh[1]) / maximum(abs.(nums_sh[1]))
@printf("   cell_A vs shipped numerator[1]:             Δrel   = %.2e   %s\n",
        relAn, relAn < 1e-10 ? "OK" : "*** MISMATCH ***")

# recover poles (roots of Q) for each cell
function poleset(Q)
    # roots of Q(z)=Σ q_k z^{k-1}; degree = length(Q)-1
    deg = length(Q) - 1
    deg == 0 && return ComplexF64[]
    # companion via Polynomials-free: use roots of reversed coeffs
    sort(abs.(roots_of(Q)))
end
# minimal root finder: companion matrix eigenvalues
using LinearAlgebra: eigvals
function roots_of(c)
    n = length(c) - 1
    n == 0 && return ComplexF64[]
    cc = ComplexF64.(c) ./ c[end]
    C = zeros(ComplexF64, n, n)
    for i = 1:n-1; C[i+1, i] = 1; end
    for i = 1:n; C[i, n] = -cc[i]; end
    eigvals(C)
end

for (nm, cl) in [("A (m,m)", cA), ("B ceil (m-1,m)", cBc), ("B floor", cBf)]
    rts = roots_of(cl.denominator)
    rr = sort(real.(rts))
    z = 0.5 + 0im
    e1 = abs(polyval(cl.numerators[1], z) / polyval(cl.denominator, z) - polyval(P1, z) / polyval(Qt, z))
    @printf("3. cell %-16s  ν=%d  roots≈%s  eval-err@0.5=%.2e  m_cur=%d n=%d K=%d g=%.1e\n",
            nm, length(cl.denominator) - 1, string(round.(rr, digits = 4)), e1,
            cl.m_cur, cl.n, cl.K, gap(cl))
end

# --- Check 2: d=1 oracles ---------------------------------------------------
println()
expjet(N) = ComplexF64[inv(factorial(big(k))) for k = 0:N]
m1 = 6
je = expjet(2m1 + 2)
cA1 = cell_A([je], m1)
cB1 = cell_B_ceil([je], m1)
rpA = robust_pade(je, m1, m1; method = :svd)
rpB = robust_pade(je, m1 - 1, m1; method = :svd)
dA = qdiff(cA1.denominator, rpA.b ./ rpA.b[1]) / maximum(abs.(rpA.b ./ rpA.b[1]))
dB = qdiff(cB1.denominator, rpB.b ./ rpB.b[1]) / maximum(abs.(rpB.b ./ rpB.b[1]))
@printf("2. d=1: cell_A vs robust_pade(%d,%d):    Δrel(Q)=%.2e  %s\n",
        m1, m1, dA, dA < 1e-8 ? "OK" : "*** MISMATCH ***")
@printf("   d=1: cell_B_ceil vs robust_pade(%d,%d): Δrel(Q)=%.2e  %s\n",
        m1 - 1, m1, dB, dB < 1e-8 ? "OK" : "*** MISMATCH ***")
@printf("   d=1: are A and B different? ν_A=%d ν_B=%d μ_A=%d μ_B=%d  (expect μ_A>μ_B)\n",
        length(cA1.denominator) - 1, length(cB1.denominator) - 1,
        length(cA1.numerators[1]) - 1, length(cB1.numerators[1]) - 1)

# --- Check 4: select on the rational (should be a near-tie, both exact) -----
println()
chosen, tbl = select([cA, cBc], jets)
@printf("4. select on rational {2,3}: R=%s g=%s K=%s -> chose :%s\n",
        string(round.(log10.(tbl.R .+ 1e-300), digits = 1)),
        string(round.(tbl.g, digits = 1)), string(tbl.K), tbl.chosen)
println()
println("smoke done.")
