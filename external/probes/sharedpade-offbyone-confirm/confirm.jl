# =============================================================================
# confirm.jl — adversarial empirical confirmation of the SharedPade.jl:117
# `_toeplitz_block` off-by-one (worklog 064 leading suspect C1).
#
# Contract under test (SP.1.1 / ADR-0019, src/SharedPade.jl:46-53): for d=1,
#   shared_denominator_pade([jet], m)  MUST reproduce  robust_pade(jet, m, m;
#   method=:svd)  bit-for-bit (modulo SVD sign), i.e. identical denominator Q.
#
# Decisive test (worklog 064 Pickup item 2): the contract is invisible on
# EXACT-RATIONAL jets (the true Q annihilates BOTH the (m,m) and the shifted
# (m-1,m) matching window) — and every shared_pade_test.jl oracle feeds exact
# rationals, which is why the suite is green. It should BITE on a truncated
# TRANSCENDENTAL jet. So we compare Q on both kinds of input.
#
# Read-only diagnostic: imports the package, computes, prints. No src edits.
# Run:  julia --project=. external/probes/sharedpade-offbyone-confirm/confirm.jl
#
# RESULT (2026-06-02, BigFloat 256-bit) — BUG CONFIRMED:
#   rational control 1/(1-z/2):  rel = 0.0  for m=2..5  (bit-identical ⇒ why
#                                the all-rational test suite is green).
#   exp(z):    rel = 16.7% / 10.0% / 7.1% / 5.6% / 4.5%   for m=2..6
#   log(1+z):  rel = 50%  / 24%  / 24%  / 24%  / 21%       for m=2..6
#   AND degP(shared) = degP(robust) - 1 on EVERY transcendental case — the
#   numerator collapses exactly one degree (a_m → 0), the precise mechanism
#   worklog 064 predicted. The d=1 contract (SP.1.1 / ADR-0019) is violated
#   on transcendental jets; adaptive per-step m ⇒ intermittent discontinuity.
# =============================================================================

using PadeTaylor: shared_denominator_pade, robust_pade
using Printf

setprecision(BigFloat, 256)

"Pad two coeff vectors to equal length with zeros, return the ∞-norm of their difference."
function qdiff(q1::AbstractVector, q2::AbstractVector)
    L = max(length(q1), length(q2))
    pad(v) = [i ≤ length(v) ? v[i] : zero(eltype(v)) for i in 1:L]
    return maximum(abs.(pad(q1) .- pad(q2)))
end

function compare(label, jet, m)
    P  = robust_pade(jet, m, m; method = :svd)        # scalar GGT oracle
    (nums, Q) = shared_denominator_pade([jet], m)     # d=1 shared-denominator
    Qr = P.b ./ P.b[1]                                # both already b[1]=1, but be safe
    Qs = Q  ./ Q[1]
    d  = qdiff(Qr, Qs)
    rel = d / maximum(abs.(Qr))
    @printf("  [%-18s m=%d]  ‖ΔQ‖∞ = %.3e   rel = %.3e   degQ: robust=%d shared=%d   degP: robust=%d shared=%d\n",
            label, m, Float64(d), Float64(rel),
            length(Qr) - 1, length(Qs) - 1, length(P.a) - 1, length(nums[1]) - 1)
    return rel
end

# --- jets -------------------------------------------------------------------
# exp(z): c_k = 1/k!  (transcendental — should BITE if the bug is real)
expjet(N)  = BigFloat[ inv(factorial(big(k))) for k in 0:N ]
# 1/(1 - z/2): c_k = (1/2)^k  (EXACT rational, pole at z=2 — control, should AGREE)
ratjet(N)  = BigFloat[ (big(1)/2)^k for k in 0:N ]
# log(1+z): c_k = (-1)^(k+1)/k for k≥1, c_0=0  (transcendental — should BITE)
logjet(N)  = BigFloat[ k == 0 ? big(0.0) : (-1)^(k+1) / big(k) for k in 0:N ]

println("=" ^ 78)
println("SharedPade d=1 vs robust_pade(:,m,m;:svd)  —  contract SP.1.1 / ADR-0019")
println("BigFloat 256-bit.  rel ≈ 0 ⇒ contract holds;  rel = O(1) ⇒ bug confirmed.")
println("=" ^ 78)

println("\nCONTROL — exact rational 1/(1 - z/2)  [expect rel ≈ 0]:")
maxrat = 0.0
for m in 2:5
    global maxrat = max(maxrat, compare("rational", ratjet(2m + 2), m))
end

println("\nTRANSCENDENTAL — exp(z)  [expect rel = O(1) iff bug real]:")
maxexp = 0.0
for m in 2:6
    global maxexp = max(maxexp, compare("exp", expjet(2m + 2), m))
end

println("\nTRANSCENDENTAL — log(1+z)  [expect rel = O(1) iff bug real]:")
maxlog = 0.0
for m in 2:6
    global maxlog = max(maxlog, compare("log1p", logjet(2m + 2), m))
end

println("\n" * "=" ^ 78)
@printf("VERDICT:  rational(control) max rel = %.2e   |   exp max rel = %.2e   log1p max rel = %.2e\n",
        maxrat, maxexp, maxlog)
if maxexp > 1e-6 || maxlog > 1e-6
    println("⇒ CONFIRMED: d=1 shared-denominator Q DIVERGES from the scalar GGT oracle")
    println("  on transcendental jets while AGREEING on exact rationals — exactly the")
    println("  worklog-064 off-by-one signature (SharedPade.jl:117).")
else
    println("⇒ NOT confirmed: denominators agree on transcendental jets too.")
end
println("=" ^ 78)
