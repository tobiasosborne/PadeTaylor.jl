# =============================================================================
# verify.jl — cross-check the block-Toeplitz-determinant oracle (oracle.jl)
# against (a) the exactly-known Q used to construct each test input, and
# (b) PadeTaylor.SharedPade.shared_denominator_pade (the SVD+QR path V1c is
# built to validate).
#
# This is the V1c analogue of external/probes/calgo766-oracle/capture.jl: an
# INDEPENDENT route's verification record.  verify.jl computes the oracle
# output, the SharedPade output, and the ground truth, and reports the
# observed agreement tolerances.  It does NOT add asserts to the project
# test suite — bead V1e wires oracle.jl into test/ separately.
#
# Inputs are the SP.1.1 (d=1 scalar rational) and SP.1.2 (d=2 shared-Q) cases
# from test/shared_pade_test.jl, fed in two ways:
#   - Float64  jets — to compare against SharedPade (a Float64 SVD path);
#   - Rational{BigInt} jets — to obtain the EXACT oracle output and compare
#     it bit-for-bit against the exactly-known Q (the "tight oracle" claim).
#
# Run (one Julia process — CLAUDE.md Rule 7):
#   julia --project=../../.. external/probes/shared-pade-determinant-oracle/verify.jl
# =============================================================================

using LinearAlgebra
using Printf
using PadeTaylor: shared_denominator_pade

include(joinpath(@__DIR__, "oracle.jl"))
using .SharedPadeDeterminantOracle: shared_q_via_determinant

# --- helpers (mirrors of test/shared_pade_test.jl) ---------------------------

# Formal power-series long division: Taylor coefficients [c_0 … c_N] of P/Q.
# Generic in element type T so it serves both Float64 and Rational{BigInt}.
function ratio_jet(P::Vector{T}, Q::Vector{T}, N::Int) where {T}
    c = zeros(T, N + 1)
    for k = 0:N
        s = k + 1 ≤ length(P) ? P[k + 1] : zero(T)
        for j = 1:k
            qj = j + 1 ≤ length(Q) ? Q[j + 1] : zero(T)
            s -= qj * c[k - j + 1]
        end
        c[k + 1] = s / Q[1]
    end
    return c
end

polyval(coef, z)       = sum(coef[k] * z^(k - 1) for k = 1:length(coef))
ratval(num, den, z)    = polyval(num, z) / polyval(den, z)

# Roots of a degree-≤2 polynomial (low-to-high coeffs), sorted, complex.
function roots_of(q)
    a, b, c = q[3], q[2], q[1]
    disc = sqrt(complex(b^2 - 4a * c))
    sort([(-b + disc) / (2a), (-b - disc) / (2a)], by = z -> (real(z), imag(z)))
end

# =============================================================================
# Test cases — SP.1.1 and SP.1.2 verbatim from test/shared_pade_test.jl.
# Each carries both a Float64 spec and a Rational{BigInt} spec of the SAME
# rational function (the SP test constants are all exact decimals).
# =============================================================================

# SP.1.1 — d=1: f = (1 + 2z) / (1 - 0.5z + 0.3z²), denominator degree m=2.
sp11 = (
    name = "SP.1.1 d=1 scalar rational",
    m    = 2,
    Pf   = [[1.0, 2.0]],
    Qf   = [1.0, -0.5, 0.3],
    Pr   = [Rational{BigInt}[1, 2]],
    Qr   = Rational{BigInt}[1, -1//2, 3//10],
)

# SP.1.2 — d=2: Q = 1 - 0.4z + 0.2z²; P1 = 1 + 0.3z; P2 = 2 - 0.7z + 0.5z².
sp12 = (
    name = "SP.1.2 d=2 shared denominator",
    m    = 2,
    Pf   = [[1.0, 0.3], [2.0, -0.7, 0.5]],
    Qf   = [1.0, -0.4, 0.2],
    Pr   = [Rational{BigInt}[1, 3//10], Rational{BigInt}[2, -7//10, 1//2]],
    Qr   = Rational{BigInt}[1, -2//5, 1//5],
)

# =============================================================================
# Run the verification.
# =============================================================================

function verify_case(case)
    println("=" ^ 74)
    println(case.name, "   (m = ", case.m, ", d = ", length(case.Pf), ")")
    println("-" ^ 74)
    m = case.m
    all_ok = true

    # --- (a) EXACT route: Rational{BigInt} jets -----------------------------
    # The oracle on exact input must return the EXACTLY-KNOWN Q (and the exact
    # numerators), bit-for-bit — the "tight oracle" property.
    jets_r = [ratio_jet(case.Pr[i], case.Qr, 2m) for i = 1:length(case.Pr)]
    nums_r, den_r = shared_q_via_determinant(jets_r, m)   # tol defaults to 0

    # Exact Q is normalised to Q(0)=1 already (Qr[1] == 1); compare directly.
    Qexact = case.Qr
    exact_den_ok = den_r == Qexact
    @printf("  EXACT (Rational{BigInt}) oracle vs known Q:  %s\n",
            exact_den_ok ? "den == known Q  EXACTLY" :
                           "MISMATCH  oracle=$(den_r)  known=$(Qexact)")
    all_ok &= exact_den_ok
    for i = 1:length(case.Pr)
        exact_num_ok = nums_r[i] == case.Pr[i]
        @printf("    component %d numerator: %s\n", i,
                exact_num_ok ? "P_$i == known P_$i  EXACTLY" :
                               "MISMATCH  oracle=$(nums_r[i])  known=$(case.Pr[i])")
        all_ok &= exact_num_ok
    end

    # --- (b) Float64 route: oracle vs SharedPade vs ground truth ------------
    jets_f = [ratio_jet(case.Pf[i], case.Qf, 2m) for i = 1:length(case.Pf)]
    nums_o, den_o = shared_q_via_determinant(jets_f, m)        # determinant oracle
    nums_s, den_s = shared_denominator_pade(jets_f, m)         # SharedPade SVD+QR

    # Denominator coefficient agreement (both Q(0)=1-normalised).
    den_coef_err_truth  = maximum(abs.(den_o .- case.Qf))
    den_coef_err_shared = maximum(abs.(den_o .- den_s))
    @printf("  Float64 oracle denominator:  |Δ vs truth| = %.3e   |Δ vs SharedPade| = %.3e\n",
            den_coef_err_truth, den_coef_err_shared)
    all_ok &= (den_coef_err_truth < 1e-12) && (den_coef_err_shared < 1e-12)

    # Pole locations: roots of the recovered denominator.
    if length(den_o) == 3
        r_oracle = roots_of(den_o)
        r_truth  = roots_of(case.Qf)
        r_shared = length(den_s) == 3 ? roots_of(den_s) : r_oracle
        pole_err_truth  = maximum(abs.(r_oracle .- r_truth))
        pole_err_shared = maximum(abs.(r_oracle .- r_shared))
        @printf("  Float64 oracle pole roots:   |Δ vs truth| = %.3e   |Δ vs SharedPade| = %.3e\n",
                pole_err_truth, pole_err_shared)
        all_ok &= (pole_err_truth < 1e-10) && (pole_err_shared < 1e-10)
    end

    # Function-value agreement — the load-bearing invariant (SP.1.1 docstring).
    zs = (0.1, -0.2, 0.33, 0.25 + 0.1im, -0.41im)
    valmax_truth  = 0.0
    valmax_shared = 0.0
    for i = 1:length(case.Pf)
        for z in zs
            v_oracle = ratval(nums_o[i], den_o, z)
            v_shared = ratval(nums_s[i], den_s, z)
            v_truth  = ratval(case.Pf[i], case.Qf, z)
            valmax_truth  = max(valmax_truth,  abs(v_oracle - v_truth))
            valmax_shared = max(valmax_shared, abs(v_oracle - v_shared))
        end
    end
    @printf("  Float64 oracle values:       |Δ vs truth| = %.3e   |Δ vs SharedPade| = %.3e\n",
            valmax_truth, valmax_shared)
    all_ok &= (valmax_truth < 1e-12) && (valmax_shared < 1e-12)

    println("-" ^ 74)
    println(all_ok ? "  CASE PASSED" : "  CASE FAILED — investigate (CLAUDE.md Rule 2)")
    return all_ok
end

# --- failure-mode check: jet too short must throw ----------------------------
function verify_throws()
    println("=" ^ 74)
    println("Failure-mode check: jet too short must throw (CLAUDE.md Rule 1)")
    println("-" ^ 74)
    short = [[1.0, 2.0]]                       # length 2, need m+1 = 4
    threw = false
    try
        shared_q_via_determinant(short, 3)
    catch e
        threw = true
        @printf("  threw as expected: %s\n",
                first(split(sprint(showerror, e), '\n')))
    end
    println(threw ? "  PASSED" : "  FAILED — oracle did not throw on a short jet")
    return threw
end

function main()
    println("\nBlock-Toeplitz-determinant shared-Q oracle — verification (bead V1c)\n")
    ok = true
    ok &= verify_case(sp11)
    ok &= verify_case(sp12)
    ok &= verify_throws()
    println("=" ^ 74)
    if ok
        println("VERIFICATION PASSED: the determinant oracle reproduces the exactly-")
        println("known Q (bit-for-bit on Rational input) and agrees with SharedPade")
        println("to < 1e-12 on Float64 input.  No SVD was used.")
    else
        println("VERIFICATION FAILED: see the per-case lines above.")
    end
    return ok
end

main() || exit(1)
