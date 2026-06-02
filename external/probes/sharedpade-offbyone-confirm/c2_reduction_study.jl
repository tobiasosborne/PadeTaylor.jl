# =============================================================================
# c2_reduction_study.jl — derive the C2 graceful-reduction / isolation criterion
# for the shared-denominator Pade, using the CORRECT (+2) window (worklog 066).
#
# Question: after the +2 window fix, the over-determined d>=2 stack can be
#   (i)  EXACT-NULL  (components share a genuine degree-m denominator): one σ≈0,
#        ρ=count(σ>τ)=m_cur, isolated 1-D null space -> use it.
#   (ii) RANK-DEFICIENT (entire / far-from-pole): many σ below τ at large m ->
#        reduce m (GGT md:136 reduce-n-to-ρ).
#   (iii) FULL COLUMN RANK (ρ=m_cur+1, NO σ≤τ): no exact null; the smallest-σ
#        right vector is a least-squares "best shared Q". Accept (if isolated &
#        Q(0) normalizable) or reduce further? -- THIS is the open criterion.
#
# We replicate the block build (+2), the reduction loop (parameterized rule),
# the QR-reweighting, and numerator recovery, and study TWO regimes:
#   A) meromorphic d=2 with KNOWN shared poles {2,3}  -> expect exact-null at m=2
#   B) entire harmonic d=2 [cos,-sin]                 -> rank-deficient, no poles
# For each rule we print the reduction trace and the resulting approximant
# accuracy at a test point, to pin the correct accept/reduce/throw rule.
#
# Read-only: no src edits, no PadeTaylor import (uses LinearAlgebra.svd on F64).
# =============================================================================
using LinearAlgebra, Printf

# +2 (GGT diagonal (m,m)) Toeplitz block: blk[r,c] = c_{m+r-c+1}, top-left c_{m+1}
function blk_plus2(c, m)
    A = zeros(Float64, m, m + 1)
    for r = 1:m, cc = 1:(m + 1)
        idx = m + r - cc + 2
        if 1 <= idx <= length(c); A[r, cc] = c[idx]; end
    end
    A
end
# degree-m numerator upper block: up[r,c] = c_{r-c}
function upper(c, m)
    U = zeros(Float64, m + 1, m + 1)
    for cc = 1:(m + 1), r = cc:(m + 1)
        idx = r - cc + 1
        if 1 <= idx <= length(c); U[r, cc] = c[idx]; end
    end
    U
end

# formal P/Q Taylor jet to order N
function ratio_jet(P, Q, N)
    c = zeros(Float64, N + 1)
    for k = 0:N
        s = k + 1 <= length(P) ? P[k+1] : 0.0
        for j = 1:k
            qj = j + 1 <= length(Q) ? Q[j+1] : 0.0
            s -= qj * c[k-j+1]
        end
        c[k+1] = s / Q[1]
    end
    c
end
polyval(c, z) = sum(c[k] * z^(k-1) for k = 1:length(c))

# parameterized shared-Pade with reduction RULE; returns (Q, nums, trace, threw)
function shared(jets, m0; rule = :ggt_strict, tol = 1e-14)
    d = length(jets)
    cnorm = norm(reduce(vcat, jets)); τ = tol * cnorm
    m = m0; trace = String[]
    local A, S, V, b
    while true
        A = reduce(vcat, [blk_plus2(c, m) for c in jets])
        F = svd(A; full = true); S = F.S; V = F.V
        ρ = count(>(τ), S)
        gap = length(S) >= 2 ? S[end-1] / max(S[end], eps()) : Inf  # σ_{m}/σ_{m+1}
        push!(trace, @sprintf("m=%d ρ=%d/%d σ_min=%.2e gap(σ_m/σ_min)=%.2e", m, ρ, m+1, S[end], gap))
        ρ == 0 && return (nothing, nothing, trace, :throw_allzero)
        accept = false
        if rule == :current            # ρ>=m break (the buggy vector rule)
            accept = ρ >= m
        elseif rule == :ggt_strict     # accept only exact isolated null (ρ==m)
            accept = ρ == m
        elseif rule == :ls_gap         # accept exact-null OR isolated least-squares (gap large)
            accept = (ρ == m) || (ρ == m + 1 && gap > 1e3)
        elseif rule == :ggt_zlambda    # ρ==m isolated null, + z^λ cancellation on Q(0)≈0
            accept = ρ == m
        end
        if accept
            b0 = V[:, end]
            D = Diagonal(abs.(b0) .+ sqrt(eps()))
            QF = qr(adjoint(A * D)); b = D * QF.Q[:, m+1]; b ./= norm(b)
            nums = [ upper(c, m) * b for c in jets ]
            if rule == :ggt_zlambda
                # MIRROR scalar RobustPade._trim_and_normalise:461-476 — cancel
                # the common z^λ factor (leading near-zeros of b = spurious pole
                # at the expansion centre) instead of throwing Q(0)≈0.
                lam_idx = findfirst(x -> abs(x) > tol, b)
                lam_idx === nothing && return (nothing, nothing, trace, :throw_allzero_b)
                lam = lam_idx - 1
                if lam > 0
                    push!(trace, "  cancel z^$lam (Q(0)≈0 leading factor)")
                    b = b[(lam+1):end]; nums = [a[(lam+1):end] for a in nums]
                end
            else
                if abs(b[1]) <= τ && abs(b[1]) <= tol
                    return (nothing, nothing, trace, :throw_Q0)
                end
            end
            Q = b ./ b[1]
            nums = [ a ./ b[1] for a in nums ]
            return (Q, nums, trace, :ok)
        end
        # not accepted: reduce
        m = (ρ < m) ? ρ : m - 1
        m < 1 && return (nothing, nothing, trace, :throw_noreduce)
    end
end

function study(name, jets, m0, testz, truef; rescale_h = 1.0)
    js = [[rescale_h^(k-1) * c[k] for k = 1:length(c)] for c in jets]
    println("\n", "="^76, "\n$name  (request m=$m0, eval at t=$testz)")
    for rule in (:current, :ggt_zlambda)
        Q, nums, trace, status = shared(js, m0; rule = rule)
        @printf("\n  RULE %-12s -> %s\n", rule, status)
        for t in trace; println("      ", t); end
        if status == :ok
            qz = polyval(Q, testz)
            errs = [abs(polyval(nums[i], testz) / qz - truef[i]) for i in 1:length(jets)]
            @printf("      settled m=%d, Q(0)=%.3f, eval err per comp: %s\n",
                    length(Q)-1, Q[1], join([@sprintf("%.2e", e) for e in errs], ", "))
        end
    end
end

# --- Regime A: meromorphic d=2, shared poles {2,3} ---------------------------
Qt = [1.0, -5/6, 1/6]                 # (1-z/2)(1-z/3)
P1 = [1.0, 2.0]; P2 = [3.0, -1.0]
m0 = 5                                  # over-ask (true shared degree is 2)
j1 = ratio_jet(P1, Qt, 2m0); j2 = ratio_jet(P2, Qt, 2m0)
tz = 0.4
study("A: meromorphic shared poles {2,3}, true Q deg 2", [j1, j2], m0, tz,
      [polyval(P1,tz)/polyval(Qt,tz), polyval(P2,tz)/polyval(Qt,tz)])

# --- Regime B: entire harmonic [cos,-sin], stepper-like large m + rescale ----
N = 30
cosj = [iseven(k) ? Float64((-1)^(k÷2)//factorial(big(k))) : 0.0 for k = 0:N]
nsinj = [isodd(k) ? Float64(-(-1)^((k-1)÷2)//factorial(big(k))) : 0.0 for k = 0:N]
h = 0.7
study("B: entire harmonic [cos,-sin], rescaled by h=0.7", [cosj, nsinj], 15, 1.0,
      [cos(h), -sin(h)]; rescale_h = h)
