# external/probes/out-of-class-guard/probe2.jl
#
# STRESS-TEST extension of probe.jl: push pole-bridge marches as close to the
# pole as the corpus legitimately does, and measure δ for an OFF-AXIS complex
# pole step (the hardest in-class case), to bound max in-class δ tightly.
# Also reproduces the EXACT solve_pade clamped-step trajectory for e^{1/z}.

using PadeTaylor
import PadeTaylor: taylor_coefficients_2nd
import PadeTaylor.RobustPade: robust_pade
import PadeTaylor.PadeStepper: PadeStepperState, pade_step_with_pade!,
                               _rescale_by_powers, _evaluate_pade

const TINY = 1e-300

function delta_at(state, f, order, h)
    coefs   = taylor_coefficients_2nd(f, state.z, state.u, state.up, order)
    scaled  = _rescale_by_powers(coefs, eltype(coefs)(h))
    m = n   = order ÷ 2
    P    = robust_pade(scaled, m, n)
    one_T = one(eltype(scaled))
    P1   = _evaluate_pade(P, one_T)
    scaled_lo = scaled[1:(end - 2)]
    P_lo  = robust_pade(scaled_lo, m - 1, n - 1)
    P1_lo = _evaluate_pade(P_lo, one_T)
    return abs(P1 - P1_lo) / (abs(P1) + TINY)
end

function march(label, f, T, u0, up0, z0, zend, h; order = 30, verbose = true)
    verbose && println("\n=== $label ===")
    state = PadeStepperState{T}(z0, u0, up0)
    # Unit direction along the (possibly complex) span; remaining distance is
    # the projection of (zend - z) onto that direction.
    udir = (zend - z0) / abs(zend - z0)
    deltas = Float64[]
    step = 0
    remaining() = real(conj(udir) * (zend - state.z))
    while remaining() > 1e-12
        step += 1
        step > 2000 && break
        h_step = udir * min(abs(h), remaining())
        local δ
        try
            δ = delta_at(state, f, order, h_step)
        catch e
            verbose && println("  step $step z=$(round(real(state.z);digits=4)) δ=THROW")
            break
        end
        push!(deltas, δ)
        try
            pade_step_with_pade!(state, f, order, h_step)
        catch e
            verbose && println("  step $step δ=$(round(δ;sigdigits=4)) ADVANCE-THROW")
            break
        end
        if verbose
            zr = T <: Complex ? round(state.z; digits=4) : round(state.z; digits=5)
            println("  step $step z=$zr δ=$(δ==0 ? "0" : round(δ;sigdigits=4))")
        end
    end
    !isempty(deltas) && verbose &&
        println("  max δ=$(round(maximum(deltas);sigdigits=4))")
    return deltas
end

allmax = Float64[]

# ℘ pole at z=1: march RIGHT up to 0.98, 0.99 at small h (closest legit approach)
fW(z, u, up) = 6 * u^2
push!(allmax, maximum(march("℘ march 0→0.99 h=0.05", fW, Float64,
    1.07182251641691742068161429481250292922453268,
    1.71033735317678620846433189149624950839231872, 0.0, 0.99, 0.05)))
push!(allmax, maximum(march("℘ march 0→0.98 h=0.02", fW, Float64,
    1.07182251641691742068161429481250292922453268,
    1.71033735317678620846433189149624950839231872, 0.0, 0.98, 0.02)))

# tan pole π/2≈1.5708: march to 1.56 at small h
ftan(z, u, up) = 2 * u + 2 * u^3
push!(allmax, maximum(march("tan march 0→1.56 h=0.05", ftan, Float64,
    0.0, 1.0, 0.0, 1.56, 0.05)))

# OFF-AXIS complex pole — tanh, pole at i·π/2 ≈ 1.5708i.  u''=-2u+2u^3, u=tanh(z).
# March in the complex direction from 0 toward i·1.5 at h=0.1i (approaches the
# imaginary-axis pole).  This is the CPB.6 family (path_network territory).
ftanh(z, u, up) = -2 * u + 2 * u^3
d_tanh = march("tanh complex pole iπ/2: march 0→1.5i h=0.1i", ftanh,
    ComplexF64, 0.0 + 0im, 1.0 + 0im, 0.0 + 0im, 0.0 + 1.5im, 0.1im)
isempty(d_tanh) || push!(allmax, maximum(d_tanh))

# sn complex pole — jacobi sn, pole at i·K'.  u''=-(5/4)u+(1/2)u^3 (CPB.7).
fsn(z, u, up) = -(5/4) * u + (1/2) * u^3
d_sn = march("sn complex: march 0→1.0i h=0.1i", fsn, ComplexF64,
    0.0 + 0im, 1.0 + 0im, 0.0 + 0im, 0.0 + 1.0im, 0.1im)
isempty(d_sn) || push!(allmax, maximum(d_sn))

println("\n################ STRESS SUMMARY ################")
println("MAX IN-CLASS δ (stress, all steps) = ", round(maximum(allmax); sigdigits = 4))

# ---------------------------------------------------------------------------
# EXACT solve_pade clamped-step trajectory for e^{1/z} (z0=-1 → -0.02, h=0.1).
# Reproduce the driver's gap-clamping to show the precise δ sequence the GUARD
# will see in the real driver.
# ---------------------------------------------------------------------------
fess(z, u, up) = u * (1 + 2z) / z^4
println("\n=== e^{1/z} EXACT driver trajectory (clamped) ===")
state = PadeStepperState{Float64}(-1.0, exp(-1.0), -exp(-1.0))
zend = -0.02; h_max = 0.1; dir = sign(zend - (-1.0))
dseq = Float64[]
while dir * (zend - state.z) > 0
    gap = zend - state.z
    h_step = dir * min(h_max, abs(gap))
    δ = delta_at(state, fess, 30, h_step)
    push!(dseq, δ)
    pade_step_with_pade!(state, fess, 30, h_step)
    println("  z=$(round(state.z;digits=4))  δ=$(δ==0 ? "0" : round(δ;sigdigits=4))")
    dir * (zend - state.z) > 0 || break
end
println("  δ sequence: ", [round(x;sigdigits=3) for x in dseq])
println("###################################################")
