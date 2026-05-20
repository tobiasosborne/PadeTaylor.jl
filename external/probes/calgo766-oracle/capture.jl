# =============================================================================
# capture.jl — ccall ACM Calgo 766 VECTOR_PADE (via the libcalgo766.so shim),
# feed it the SharedPade test inputs, pin the type-II outputs to oracles.txt,
# and sanity-check against SharedPade.shared_denominator_pade.
#
# This is the v0.2 analogue of external/probes/padeapprox-oracle/capture.m:
# an INDEPENDENT reference implementation (Cabay-Jones-Labahn 1997, ACM TOMS
# Algorithm 766) against which bead V1e will mutation-prove our own shared-Q
# code.  capture.jl computes nothing itself — it marshals inputs to FORTRAN
# and marshals the answer back.
#
# Run (one Julia process — CLAUDE.md Rule 7):
#   julia --project=../../.. external/probes/calgo766-oracle/capture.jl
# (build.sh must have produced libcalgo766.so first.)
#
# CONVENTION MAPPING (see README.md and wrapper.f90 header for full detail):
#   Calgo type-II row-0 of S_star gives, with A(0):=1 and A(beta):=f_beta,
#       Q      = S_star(:,0,0)        f_beta ~= P_beta / Q
#       P_beta = S_star(:,0,beta)
#   coefficients low-to-high, NOT normalised to Q(0)=1 — we renormalise here
#   exactly as SharedPade does (divide P and Q by Q[1]).
# =============================================================================

using LinearAlgebra
using Printf
using PadeTaylor: shared_denominator_pade

const LIB = joinpath(@__DIR__, "libcalgo766.so")
isfile(LIB) || error("capture.jl: $LIB missing — run ./build.sh first.")

# --- ccall shim --------------------------------------------------------------
# vpade_type2(k, nvec, nrm, Amat, tau, Qout, Pout, flag) — all by reference,
# arrays column-major.  Amat is (nrm+1)x(k+1); Pout is (nrm+1)xk.
function calgo766_type2(series::Vector{Vector{Float64}}, nvec::Vector{Int};
                        tau::Float64 = 1.0e3)
    k    = length(series) - 1                 # k+1 series; k = #components d
    nrm  = sum(nvec)
    @assert length(nvec) == k + 1
    @assert all(length(s) ≥ nrm + 1 for s in series) "series too short for ||n||"

    Amat = zeros(Float64, nrm + 1, k + 1)
    for j in 1:(k + 1), l in 1:(nrm + 1)
        Amat[l, j] = series[j][l]
    end
    Qout = zeros(Float64, nrm + 1)
    Pout = zeros(Float64, nrm + 1, k)
    flag = Ref{Cint}(0)

    ccall((:vpade_type2, LIB), Cvoid,
          (Ref{Cint}, Ptr{Cint}, Ref{Cint}, Ptr{Float64}, Ref{Float64},
           Ptr{Float64}, Ptr{Float64}, Ref{Cint}),
          Cint(k), Cint.(nvec), Cint(nrm), Amat, tau, Qout, Pout, flag)

    return (Q = Qout, P = Pout, flag = Int(flag[]))
end

# Trim trailing near-zeros and renormalise so Q(0)=1, matching SharedPade.
function normalise(Q::Vector{Float64}, P::Matrix{Float64}; trimtol = 1e-9)
    q0 = Q[1]
    abs(q0) > trimtol || error("capture.jl: Calgo Q(0)=$q0 ~ 0; cannot renormalise.")
    Qn = Q ./ q0
    lastq = findlast(x -> abs(x) > trimtol, Qn)
    Qn = Qn[1:lastq]
    d  = size(P, 2)
    Pn = Vector{Vector{Float64}}(undef, d)
    for j in 1:d
        col   = P[:, j] ./ q0
        lastp = findlast(x -> abs(x) > trimtol, col)
        Pn[j] = lastp === nothing ? [0.0] : col[1:lastp]
    end
    return (Pn, Qn)
end

# Formal power-series long division — the same helper test/shared_pade_test.jl
# uses, so capture.jl builds its jets the identical way.
function ratio_jet(P::Vector{Float64}, Q::Vector{Float64}, N::Int)
    c = zeros(Float64, N + 1)
    for k in 0:N
        s = k + 1 ≤ length(P) ? P[k + 1] : 0.0
        for j in 1:k
            qj = j + 1 ≤ length(Q) ? Q[j + 1] : 0.0
            s -= qj * c[k - j + 1]
        end
        c[k + 1] = s / Q[1]
    end
    return c
end

polyval(c, z)         = sum(c[k] * z^(k - 1) for k in 1:length(c))
ratval(num, den, z)   = polyval(num, z) / polyval(den, z)

# --- write Julia-readable literals ------------------------------------------
function write_vec(io, name, v)
    print(io, name, " = [")
    for (k, x) in enumerate(v)
        k > 1 && print(io, ", ")
        @printf(io, "%.18e", x)
    end
    println(io, "]")
end

# =============================================================================
# Test cases — the SharedPade SP.1.1 / SP.1.2 inputs.
#
# DEGREE CHOICE.  Calgo type-II row-0 degrees are ||n||-n(0) for Q and
# ||n||-n(beta) for P_beta.  We pick n so those degrees comfortably exceed the
# true rational's degrees: an exact rational is recovered exactly by any
# Pade of high enough type, so the oracle output is the true (P_i, Q) regardless
# of the exact n — only the *jet length* must reach ||n||+1.  We pin n in
# oracles.txt so the capture is fully reproducible.
# =============================================================================

cases = []

# --- SP.1.1 — d=1 scalar rational jet ---------------------------------------
# f = (1 + 2z) / (1 - 0.5z + 0.3z^2), type (1,2).
let
    Pnum = [1.0, 2.0]
    Qden = [1.0, -0.5, 0.3]
    k    = 1                          # 1 component  (k+1 = 2 series)
    nvec = [1, 2]                     # ||n||=3; Q deg = 3-1 = 2, P deg = 3-2 = 1
    nrm  = sum(nvec)
    A0   = [l == 1 ? 1.0 : 0.0 for l in 1:(nrm + 1)]   # reference series == 1
    f1   = ratio_jet(Pnum, Qden, nrm)
    push!(cases, (name = "SP_1_1_d1_scalar", Pnum = [Pnum], Qden = Qden,
                  series = [A0, f1], nvec = nvec, k = k))
end

# --- SP.1.2 — d=2 genuine shared denominator --------------------------------
# Q = 1 - 0.4z + 0.2z^2; P1 = 1 + 0.3z; P2 = 2 - 0.7z + 0.5z^2.
let
    Qden = [1.0, -0.4, 0.2]
    P1   = [1.0, 0.3]
    P2   = [2.0, -0.7, 0.5]
    k    = 2                          # 2 components (k+1 = 3 series)
    nvec = [2, 3, 2]                  # ||n||=7; Q deg = 5, P1 deg = 4, P2 deg = 5
    nrm  = sum(nvec)
    A0   = [l == 1 ? 1.0 : 0.0 for l in 1:(nrm + 1)]
    f1   = ratio_jet(P1, Qden, nrm)
    f2   = ratio_jet(P2, Qden, nrm)
    push!(cases, (name = "SP_1_2_d2_shared_Q", Pnum = [P1, P2], Qden = Qden,
                  series = [A0, f1, f2], nvec = nvec, k = k))
end

# =============================================================================
# Run Calgo 766, pin oracles, sanity-check against SharedPade.
# =============================================================================

function run_capture(cases)
io = open(joinpath(@__DIR__, "oracles.txt"), "w")
println(io, "# Calgo-766 type-II (simultaneous-Pade / shared-Q) oracle outputs.")
println(io, "# Source: ACM TOMS Algorithm 766 (Cabay-Jones-Labahn 1997),")
println(io, "#         DOI 10.1145/244768.244790, via libcalgo766.so.")
println(io, "# Calgo Sp/ build uses single-precision `real`: accuracy floor ~1e-6.")
println(io, "# Q, P_i coefficients low-to-high, renormalised so Q(0)=1.")
println(io, "# Format: Julia-readable literals.\n")

println("=" ^ 72)
all_ok = true

for c in cases
    res = calgo766_type2(c.series, c.nvec)
    @printf("%-22s  Calgo flag = %d\n", c.name, res.flag)
    # flag 1 = ill-conditioned (kappa >= tau) is acceptable here: the row-0
    # type-II approximant is still valid per the Calgo header.  flag 3 = bad
    # input, flag 2 = singular Sylvester matrix would be a real problem.
    res.flag == 3 && error("capture.jl: Calgo flag=3 (bad input) for $(c.name).")

    Pn, Qn = normalise(res.Q, res.P)

    # --- pin to oracles.txt --------------------------------------------------
    println(io, "# Case $(c.name) — Calgo n = $(c.nvec), k = $(c.k), flag = $(res.flag)")
    write_vec(io, "$(c.name)_Q", Qn)
    for (j, p) in enumerate(Pn)
        write_vec(io, "$(c.name)_P$(j)", p)
    end
    println(io, "$(c.name)_flag = $(res.flag)\n")

    # --- sanity check vs SharedPade -----------------------------------------
    # SharedPade uses one shared degree m; the SP tests use m=2.  We compare
    # the recovered RATIONAL FUNCTIONS (Calgo vs SharedPade vs ground truth)
    # at sample points — the same value-invariant the SP tests assert, robust
    # to the differing internal degree conventions and to sign/scale.
    m = 2
    jets = c.series[2:end]               # drop the reference series A(0)=1
    sp_nums, sp_den = shared_denominator_pade(jets, m)

    zs = (0.1, -0.2, 0.27, 0.15 + 0.1im)
    casemax = 0.0
    for (j, Ptrue) in enumerate(c.Pnum)
        for z in zs
            v_calgo  = ratval(Pn[j], Qn, z)
            v_shared = ratval(sp_nums[j], sp_den, z)
            v_true   = ratval(Ptrue, c.Qden, z)
            casemax  = max(casemax, abs(v_calgo - v_shared), abs(v_calgo - v_true))
        end
    end
    # Pole set: roots of the shared denominator (degree-2 true Q).
    rootcmp = ""
    if length(Qn) ≥ 3 && length(sp_den) == 3
        rts(q) = begin
            a, b, cc = q[3], q[2], q[1]
            disc = sqrt(complex(b^2 - 4a * cc))
            sort([(-b + disc) / (2a), (-b - disc) / (2a)], by = z -> (real(z), imag(z)))
        end
        # Calgo Q here is degree 5 (exact rational embedded high) — its two
        # genuine roots match Q's; we instead compare SharedPade roots to truth
        # and report the Calgo value-agreement above as the load-bearing check.
        rooterr = maximum(abs.(rts(sp_den) .- rts(c.Qden)))
        rootcmp = @sprintf(" | SharedPade pole-root err vs truth = %.2e", rooterr)
    end

    ok = casemax < 1e-5
    all_ok &= ok
    @printf("  max |Calgo - SharedPade| & |Calgo - truth| over sample pts = %.3e  [%s]%s\n",
            casemax, ok ? "AGREE" : "DISAGREE", rootcmp)
end

close(io)
println("=" ^ 72)
println(all_ok ?
    "SANITY CHECK PASSED: Calgo 766 agrees with SharedPade to < 1e-5 (single-precision floor)." :
    "SANITY CHECK FAILED: investigate convention mismatch vs genuine bug (CLAUDE.md Rule 2).")
println("Wrote ", joinpath(@__DIR__, "oracles.txt"))
return all_ok
end

run_capture(cases) || exit(1)
